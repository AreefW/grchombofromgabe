/*
 * Rebox level 0 of a GRChombo/Chombo checkpoint.
 *
 * This version preserves the real checkpoint structure written by
 * Chombo::AMR::writeCheckpointFile plus GRChombo's GRAMRLevel checkpoint
 * methods:
 *
 *   / attributes:
 *       max_level, num_levels, iteration, time,
 *       regrid_interval_*, steps_since_regrid_*,
 *       num_components, component_*
 *   /level_N attributes:
 *       ref_ratio, tag_buffer_size, dx, dt, time, prob_domain,
 *       is_periodic_*
 *   /level_N boxes and data
 *
 * Only level_0's boxes/data are changed. All root and level headers are copied
 * from the input checkpoint. Finer-level layouts/data are copied through
 * unchanged.
 *
 * Usage:
 *   rebox_checkpoint_level0 in_chk.hdf5 out_chk.hdf5 max_box_size [block_factor]
 */

#include <cstdlib>
#include <string>

#ifdef CH_MPI
#include "mpi.h"
#endif

#ifndef CH_USE_HDF5
#error "rebox_checkpoint_level0 requires Chombo built with CH_USE_HDF5"
#endif

#include "BRMeshRefine.H"
#include "CH_HDF5.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "IntVect.H"
#include "Interval.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "MayDay.H"
#include "ProblemDomain.H"
#include "Vector.H"
#include "parstream.H"
#include "UsingNamespace.H"

namespace
{
void usage(const char *prog)
{
    pout() << "Usage: " << prog
           << " in_chk.hdf5 out_chk.hdf5 max_box_size [block_factor]\n";
}

int parsePositiveInt(const char *text, const char *name)
{
    char *end = NULL;
    const long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value <= 0)
    {
        pout() << "Invalid " << name << ": " << text << "\n";
        MayDay::Error("invalid command-line argument");
    }
    return static_cast<int>(value);
}

string levelGroupName(const int level)
{
    char level_str[40];
    sprintf(level_str, "/level_%d", level);
    return string(level_str);
}

ProblemDomain problemDomainFromHeader(const HDF5HeaderData &header)
{
    if (header.m_box.find("prob_domain") == header.m_box.end())
    {
        MayDay::Error("level header does not contain prob_domain");
    }

    bool is_periodic[SpaceDim] = {false};
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
        char dir_str[20];
        sprintf(dir_str, "%d", dir);
        const string periodic_label = string("is_periodic_") + dir_str;
        if (header.m_int.find(periodic_label) != header.m_int.end())
        {
            is_periodic[dir] = (header.m_int.find(periodic_label)->second != 0);
        }
    }

    return ProblemDomain(header.m_box.find("prob_domain")->second,
                         is_periodic);
}

DisjointBoxLayout makeLayout(const Vector<Box> &boxes,
                             const ProblemDomain &domain)
{
    Vector<int> proc_assign;
    const int status = LoadBalance(proc_assign, boxes);
    if (status != 0)
    {
        MayDay::Error("LoadBalance failed");
    }

    DisjointBoxLayout grids(boxes, proc_assign, domain);
    grids.close();
    return grids;
}

DisjointBoxLayout makeSplitLayout(const ProblemDomain &domain,
                                  const int max_box_size,
                                  const int block_factor)
{
    Vector<Box> new_boxes;
    domainSplit(domain, new_boxes, max_box_size, block_factor);

    if (new_boxes.size() == 0)
    {
        MayDay::Error("domainSplit produced no level-0 boxes");
    }

    return makeLayout(new_boxes, domain);
}

HDF5HeaderData readDataAttributes(HDF5Handle &handle)
{
    const string group = handle.getGroup();
    HDF5HeaderData data_header;

    if (handle.setGroup(group + "/data_attributes") != 0)
    {
        MayDay::Error("checkpoint level does not contain data_attributes");
    }

    data_header.readFromFile(handle);
    handle.setGroup(group);
    return data_header;
}

// Ghost cells stored in the file around each box.
IntVect fileGhost(HDF5HeaderData &data_header)
{
    if (data_header.m_intvect.find("outputGhost") !=
        data_header.m_intvect.end())
    {
        return data_header.m_intvect["outputGhost"];
    }

    if (data_header.m_intvect.find("ghost") != data_header.m_intvect.end())
    {
        return data_header.m_intvect["ghost"];
    }

    return IntVect::Zero;
}

// Ghost cells of the LevelData in memory (what Chombo's read would allocate).
IntVect memoryGhost(HDF5HeaderData &data_header)
{
    if (data_header.m_intvect.find("ghost") != data_header.m_intvect.end())
    {
        return data_header.m_intvect["ghost"];
    }

    return fileGhost(data_header);
}

// Maximum number of Reals per H5Dread call. HDF5 1.12's MPI-IO driver passes
// the byte count to MPI as an int, so a single read must stay below 2 GiB.
const long long max_read_reals = 1LL << 26;

// Read a contiguous range of the flattened data array in pieces small enough
// for the MPI-IO driver.
void readRangeChunked(hid_t dataset, hid_t filespace, const long long start,
                      const long long count, Real *buffer)
{
    long long done = 0;
    while (done < count)
    {
        const long long piece =
            (count - done < max_read_reals) ? (count - done) : max_read_reals;

        hsize_t file_start = static_cast<hsize_t>(start + done);
        hsize_t file_count = static_cast<hsize_t>(piece);
        if (H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &file_start, NULL,
                                &file_count, NULL) < 0)
        {
            MayDay::Error("failed to select hyperslab in level data");
        }

        hid_t memspace = H5Screate_simple(1, &file_count, NULL);
        if (H5Dread(dataset, H5T_NATIVE_REAL, memspace, filespace,
                    H5P_DEFAULT, buffer + done) < 0)
        {
            MayDay::Error("failed to read level data chunk");
        }
        H5Sclose(memspace);

        done += piece;
    }
}

// Replacement for Chombo's read<FArrayBox>(..., "data", ...). Chombo reads
// each box with one H5Dread, which fails for boxes larger than 2 GiB (e.g. a
// single 256^3 level-0 box). This reads each box in chunks instead.
LevelData<FArrayBox> *readLevelData(HDF5Handle &handle,
                                    const DisjointBoxLayout &grids,
                                    const Vector<Box> &file_boxes)
{
    HDF5HeaderData data_header = readDataAttributes(handle);
    if (data_header.m_int.find("comps") == data_header.m_int.end())
    {
        MayDay::Error("data_attributes does not contain comps");
    }
    const int num_comps = data_header.m_int["comps"];
    const IntVect file_ghost = fileGhost(data_header);
    const IntVect mem_ghost = memoryGhost(data_header);

    LevelData<FArrayBox> *data =
        new LevelData<FArrayBox>(grids, num_comps, mem_ghost);

    hid_t offsets_set = H5Dopen2(handle.groupID(), "data:offsets=0",
                                 H5P_DEFAULT);
    if (offsets_set < 0)
    {
        MayDay::Error("failed to open data:offsets=0");
    }
    Vector<long long> offsets(file_boxes.size() + 1);
    if (H5Dread(offsets_set, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                &offsets[0]) < 0)
    {
        MayDay::Error("failed to read data:offsets=0");
    }
    H5Dclose(offsets_set);

    hid_t dataset = H5Dopen2(handle.groupID(), "data:datatype=0",
                             H5P_DEFAULT);
    if (dataset < 0)
    {
        MayDay::Error("failed to open data:datatype=0");
    }
    hid_t filespace = H5Dget_space(dataset);

    for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
    {
        const Box &box = grids[dit];

        int file_index = -1;
        for (int i = 0; i < file_boxes.size(); ++i)
        {
            if (file_boxes[i] == box)
            {
                file_index = i;
                break;
            }
        }
        if (file_index < 0)
        {
            MayDay::Error("layout box not found in checkpoint boxes");
        }

        const Box file_box = grow(box, file_ghost);
        const long long start = offsets[file_index];
        const long long count = offsets[file_index + 1] - start;
        if (count != static_cast<long long>(file_box.numPts()) * num_comps)
        {
            MayDay::Error("level data size does not match box, comps and "
                          "ghost cells");
        }

        FArrayBox &fab = (*data)[dit];
        if (fab.box() == file_box)
        {
            // FArrayBox memory has the same component-major ordering as the
            // file, so read straight into it.
            readRangeChunked(dataset, filespace, start, count, fab.dataPtr(0));
        }
        else
        {
            FArrayBox file_fab(file_box, num_comps);
            readRangeChunked(dataset, filespace, start, count,
                             file_fab.dataPtr(0));
            fab.setVal(0.0);
            const Box overlap = fab.box() & file_box;
            fab.copy(file_fab, overlap, 0, overlap, 0, num_comps);
        }
    }

    H5Sclose(filespace);
    H5Dclose(dataset);
    return data;
}

// Chombo writes each box with a single H5Dwrite, so boxes larger than 2 GiB
// will fail on write with the same MPI-IO limit.
void warnIfBoxesTooLargeToWrite(const LevelData<FArrayBox> &data,
                                const IntVect &output_ghost, const int level)
{
    const DisjointBoxLayout &grids = data.disjointBoxLayout();
    for (LayoutIterator lit = grids.layoutIterator(); lit.ok(); ++lit)
    {
        const long long reals =
            static_cast<long long>(grow(grids[lit], output_ghost).numPts()) *
            data.nComp();
        if (reals * static_cast<long long>(sizeof(Real)) > 2147483647LL)
        {
            pout() << "Warning: level " << level << " has box " << grids[lit]
                   << " larger than 2 GiB; writing it may fail with HDF5 "
                      "MPI-IO. Use a smaller max_box_size.\n";
            return;
        }
    }
}

LevelData<FArrayBox> *copyToNewLayout(const LevelData<FArrayBox> &old_data,
                                      const DisjointBoxLayout &new_grids)
{
    LevelData<FArrayBox> *new_data =
        new LevelData<FArrayBox>(new_grids, old_data.nComp(),
                                 old_data.ghostVect());

    const Interval comps(0, old_data.nComp() - 1);
    old_data.copyTo(comps, *new_data, comps);
    new_data->exchange();
    return new_data;
}

void readCheckpointLevel(HDF5Handle &in_handle, const int level,
                         HDF5HeaderData &level_header,
                         Vector<Box> &boxes,
                         LevelData<FArrayBox> *&data,
                         IntVect &output_ghost)
{
    in_handle.setGroup(levelGroupName(level));

    if (level_header.readFromFile(in_handle) != 0)
    {
        MayDay::Error("failed to read level header");
    }

    if (read(in_handle, boxes) != 0)
    {
        MayDay::Error("failed to read level boxes");
    }

    const ProblemDomain domain = problemDomainFromHeader(level_header);
    const DisjointBoxLayout grids = makeLayout(boxes, domain);

    HDF5HeaderData data_header = readDataAttributes(in_handle);
    output_ghost = fileGhost(data_header);
    data = readLevelData(in_handle, grids, boxes);
}

void writeCheckpointLevel(HDF5Handle &out_handle, const int level,
                          HDF5HeaderData &level_header,
                          const LevelData<FArrayBox> &data,
                          const IntVect &output_ghost)
{
    out_handle.setGroup(levelGroupName(level));

    if (level_header.writeToFile(out_handle) != 0)
    {
        MayDay::Error("failed to write level header");
    }

    if (write(out_handle, data.boxLayout()) != 0)
    {
        MayDay::Error("failed to write level boxes");
    }

    if (write(out_handle, data, "data", output_ghost) != 0)
    {
        MayDay::Error("failed to write level data");
    }
}
} // namespace

int main(int argc, char **argv)
{
#ifdef CH_MPI
    MPI_Init(&argc, &argv);
#endif

    int return_code = 0;

    if (argc < 4 || argc > 5)
    {
        usage(argv[0]);
        return_code = 2;
    }
    else
    {
        const string input_file = argv[1];
        const string output_file = argv[2];
        const int max_box_size = parsePositiveInt(argv[3], "max_box_size");
        const int block_factor =
            (argc == 5) ? parsePositiveInt(argv[4], "block_factor") : 1;

        HDF5Handle in_handle(input_file.c_str(), HDF5Handle::OPEN_RDONLY);
        HDF5Handle out_handle(output_file.c_str(), HDF5Handle::CREATE);

        HDF5HeaderData root_header;
        if (root_header.readFromFile(in_handle) != 0)
        {
            MayDay::Error("failed to read root checkpoint header");
        }

        if (root_header.m_int.find("max_level") == root_header.m_int.end())
        {
            MayDay::Error("input checkpoint root header has no max_level");
        }
        if (root_header.m_int.find("num_levels") == root_header.m_int.end())
        {
            MayDay::Error("input checkpoint root header has no num_levels");
        }

        const int num_levels = root_header.m_int["num_levels"];
        if (num_levels < 1)
        {
            MayDay::Error("input checkpoint has no levels");
        }

        pout() << "Copying root checkpoint header with max_level="
               << root_header.m_int["max_level"]
               << ", num_levels=" << num_levels << "\n";

        if (root_header.writeToFile(out_handle) != 0)
        {
            MayDay::Error("failed to write root checkpoint header");
        }

        for (int level = 0; level < num_levels; ++level)
        {
            HDF5HeaderData level_header;
            Vector<Box> boxes;
            LevelData<FArrayBox> *old_data = NULL;
            IntVect output_ghost = IntVect::Zero;

            readCheckpointLevel(in_handle, level, level_header, boxes,
                                old_data, output_ghost);

            const ProblemDomain domain = problemDomainFromHeader(level_header);
            LevelData<FArrayBox> *out_data = old_data;

            if (level == 0)
            {
                pout() << "Original level 0 boxes: " << boxes.size() << "\n";
                const DisjointBoxLayout out_grids =
                    makeSplitLayout(domain, max_box_size, block_factor);
                out_data = copyToNewLayout(*old_data, out_grids);
                pout() << "New level 0 boxes: " << out_grids.size() << "\n";
            }

            warnIfBoxesTooLargeToWrite(*out_data, output_ghost, level);
            writeCheckpointLevel(out_handle, level, level_header, *out_data,
                                 output_ghost);

            if (level == 0)
            {
                delete out_data;
            }
            delete old_data;
        }

        in_handle.close();
        out_handle.close();

        pout() << "Wrote " << output_file << "\n";
    }

#ifdef CH_MPI
    MPI_Finalize();
#endif

    return return_code;
}