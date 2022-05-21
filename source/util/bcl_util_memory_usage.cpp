// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "util/bcl_util_memory_usage.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <unistd.h>

#if defined(_WIN32)
  #if defined( psapi_FOUND)
    // to gain memory tracking on windows, link with the psapi library
    #include <windows.h>
    #define PSAPI_VERSION 1
    #include <psapi.h>
  #endif
#elif defined(__APPLE__)
  #include <iostream>
  #include <mach/mach.h>
  #include <stdint.h>
  #include <sys/sysctl.h>
  #include <unistd.h>
#endif

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, determines values if they are available
    MemoryUsage::MemoryUsage() :
      m_PeakVirtual( GetUndefined< size_t>()),
      m_Virtual( GetUndefined< size_t>()),
      m_PeakRAM( GetUndefined< size_t>()),
      m_RAM( GetUndefined< size_t>())
    {
      const size_t bytes_per_megabyte( 1 << 20);
  #if defined(_WIN32)
    #if defined( psapi_FOUND)
      // windows uses PROCESS_MEMORY_COUNTERS; which require linking with psapi.lib
      PROCESS_MEMORY_COUNTERS mem_counter;
      if( GetProcessMemoryInfo( GetCurrentProcess(), &mem_counter, sizeof( PROCESS_MEMORY_COUNTERS)))
      {
        m_Virtual     = mem_counter.PagefileUsage / bytes_per_megabyte;
        m_PeakVirtual = mem_counter.PeakPagefileUsage / bytes_per_megabyte;
        m_PeakRAM     = mem_counter.PeakWorkingSetSize / bytes_per_megabyte;
        m_RAM         = mem_counter.WorkingSetSize / bytes_per_megabyte;
      }
    #endif
  #elif defined(__APPLE__)
      // on apple, getrusage returns the peak memory
      rusage usage;
      getrusage( RUSAGE_SELF, &usage);
      m_PeakRAM = usage.ru_maxrss / bytes_per_megabyte;

      // task_info must be used to get the RAM and virtual memory size
      struct task_basic_info_64 task_information;
      mach_msg_type_number_t count( TASK_BASIC_INFO_64_COUNT);

      task_info
      (
        mach_task_self(),
        TASK_BASIC_INFO_64,
        task_info_t( &task_information),
        &count
      );

      m_RAM = task_information.resident_size / bytes_per_megabyte;

      // the virtual memory size on apple is obtained from task_information.virtual_size
      // note that it includes a large shared block (usually around 630MB)
      m_Virtual = task_information.virtual_size / bytes_per_megabyte;

  #else
      // on unix systems, there are several ways to get the memory consumption of a process
      // mallinfo has the downside that it does not work above 4GB and is generally deprecated
      // rusage/getrusage has the downside that Linux does not insert the memory values into the rusage struct
      // malloc_stats has the downside that it's a compiler-specific extension and can only be used to
      //   write the statistics to std::cerr, because it's format and level of information are likely to change
      // the way chosen here is to parse the same file parsed by top, htop, ps, and similar commands; which
      // live at /proc/{pid}/statm, /proc/{pid}/stat, /proc/{pid}/status.
      // /proc/{pid}/status was selected as the best choice because it has labeled fields, so new fields can be added
      // without changing the code below, and because it has the peak virtual and RAM levels
      io::IFStream input;

      const std::string filename( "/proc/" + util::Format()( getpid()) + "/status");
      // try to open the file at /proc/process-id/status.  This contains all the information desired
      if( !io::File::TryOpenIFStream( input, filename))
      {
        // could not open file, just return
        return;
      }
      // load in all the lines
      storage::Vector< std::string> lines( StringLineListFromIStream( input));
      // close the file stream
      io::File::CloseClearFStream( input);

      // set the field to value map to undefined for all the fields needed for this object
      storage::Map< std::string, size_t> field_to_value;
      field_to_value[ "VmPeak"] = util::GetUndefined< size_t>();
      field_to_value[ "VmSize"] = util::GetUndefined< size_t>();
      field_to_value[ "VmHWM"]  = util::GetUndefined< size_t>();
      field_to_value[ "VmRSS"]  = util::GetUndefined< size_t>();

      for
      (
        storage::Vector< std::string>::const_iterator itr( lines.Begin()), itr_end( lines.End());
        itr != itr_end;
        ++itr
      )
      {
        // find the position of the first ':'
        const size_t field_delimiter_pos( itr->find( ':'));

        // if the field delimiter was not found, skip this line
        if( field_delimiter_pos >= itr->size())
        {
          continue;
        }

        // find the first non-space character on the line
        const size_t field_start( itr->find_first_not_of( " \t"));

        // get the field
        const std::string field( itr->substr( field_start, field_delimiter_pos - field_start));

        // check if the field is needed
        if( !field_to_value.Has( field))
        {
          // unused field, skip to the next line
          continue;
        }

        // field was desired, parse the value
        std::string memory_value_with_unit( TrimString( itr->substr( field_delimiter_pos + 1)));

        // get the length of the memory size
        const size_t memory_size_num_chars( LengthOfUnsignedIntegerType( memory_value_with_unit));

        // load the memory size
        size_t memory_size
        (
          ConvertStringToNumericalValue< size_t>( memory_value_with_unit.substr( 0, memory_size_num_chars))
        );

        // get the unit, if there was one; otherwise, assume it is bytes
        if( memory_size_num_chars != memory_value_with_unit.size())
        {
          // look at the first letter after the number
          const size_t first_unit_letter_index( memory_value_with_unit.find_first_not_of( " \t", memory_size_num_chars));

          // cast it to lower case
          const char first_unit_letter( tolower( int( memory_value_with_unit[ first_unit_letter_index])));

          // convert the units into mb
          switch( first_unit_letter)
          {
            case 'b':
              // bytes, convert into mb
              memory_size /= bytes_per_megabyte;
              break;
            case 'k':
              // kilobytes
              memory_size /= size_t( 1024);
              break;
            case 'm':
              // megabytes
              break;
            case 'g':
              memory_size *= size_t( 1024);
              // gigabytes
              break;
            case 't':
              memory_size *= bytes_per_megabyte;
              // terabytes
              break;
            default:
              // unknown units
              BCL_MessageCrt
              (
                "Cannot read memory value: " + memory_value_with_unit
              );
              continue;
              break;
          }
        }
        else
        {
          // assume the value was in bytes, divide by bytes_per_megabyte to convert to mb
          memory_size /= bytes_per_megabyte;
        }
        field_to_value[ field] = memory_size;
      }

      // load the values for the desired fields into the appropriate size_t
      m_PeakVirtual = field_to_value[ "VmPeak"];
      m_Virtual = field_to_value[ "VmSize"];
      m_PeakRAM = field_to_value[ "VmHWM"];
      m_RAM = field_to_value[ "VmRSS"];
  #endif
    }

    //! @brief Clone function
    //! @return pointer to new MemoryUsage
    MemoryUsage *MemoryUsage::Clone() const
    {
      return new MemoryUsage( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MemoryUsage::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the peak virtual memory used, in megabytes
    //! @return the peak virtual memory used, in megabytes
    size_t MemoryUsage::GetPeakVirtualMemoryUsed() const
    {
      return m_PeakVirtual;
    }

    //! @brief get the virtual memory in use, in megabytes
    //! @return the virtual memory in use, in megabytes
    size_t MemoryUsage::GetVirtualMemoryInUse() const
    {
      return m_Virtual;
    }

    //! @brief get the peak RAM used, in megabytes
    //! @return the peak RAM used, in megabytes
    size_t MemoryUsage::GetPeakRAMUsed() const
    {
      return m_PeakRAM;
    }

    //! @brief get the RAM in use, in megabytes
    //! @return the RAM in use, in megabytes
    size_t MemoryUsage::GetRAMInUse() const
    {
      return m_RAM;
    }

    //! @brief write peak memory usage, if available, to an output stream
    //! @param STREAM stream to write the output to
    //! @return stream reference
    std::ostream &MemoryUsage::WritePeakMemoryUsageInfo( std::ostream &STREAM)
    {
      MemoryUsage usage;
      if( IsDefined( usage.m_PeakVirtual))
      {
        STREAM << "Peak virtual memory used: " << usage.m_PeakVirtual << " MB\n";
      }
      if( IsDefined( usage.m_PeakRAM))
      {
        STREAM << "Peak RAM used: " << usage.m_PeakRAM << " MB\n";
      }
      return STREAM;
    }

    //! @brief write peak memory usage, if available, to an output stream
    //! @param STREAM stream to write the output to
    //! @return stream reference
    std::ostream &MemoryUsage::WriteCurrentMemoryUsageInfo( std::ostream &STREAM)
    {
      MemoryUsage usage;
      if( IsDefined( usage.m_Virtual))
      {
        STREAM << "Virtual memory used: " << usage.m_Virtual << " MB\n";
      }
      if( IsDefined( usage.m_RAM))
      {
        STREAM << "RAM used: " << usage.m_RAM << " MB\n";
      }
      return STREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MemoryUsage::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_PeakVirtual, ISTREAM);
      io::Serialize::Read( m_Virtual, ISTREAM);
      io::Serialize::Read( m_PeakRAM, ISTREAM);
      io::Serialize::Read( m_RAM, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MemoryUsage::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_PeakVirtual, OSTREAM, INDENT);
      io::Serialize::Write( m_Virtual, OSTREAM, INDENT);
      io::Serialize::Write( m_PeakRAM, OSTREAM, INDENT);
      io::Serialize::Write( m_RAM, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace util
} // namespace bcl
