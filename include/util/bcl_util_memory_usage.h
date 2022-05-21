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

#ifndef BCL_UTIL_MEMORY_USAGE_H_
#define BCL_UTIL_MEMORY_USAGE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MemoryUsage
    //! @brief determines and holds memory usage for this process
    //!
    //! @see @link example_util_memory_usage.cpp @endlink
    //! @author mendenjl
    //! @date Sep 26, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MemoryUsage :
      public ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      size_t m_PeakVirtual; //!< peak virtual memory used, in megabytes
      size_t m_Virtual;     //!< virtual memory used currently, in megabytes
      size_t m_PeakRAM;     //!< peak ram used, in megabytes
      size_t m_RAM;         //!< ram currently used

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, determines values if they are available
      MemoryUsage();

      //! @brief Clone function
      //! @return pointer to new MemoryUsage
      MemoryUsage *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the peak virtual memory used, in megabytes
      //! @return the peak virtual memory used, in megabytes
      size_t GetPeakVirtualMemoryUsed() const;

      //! @brief get the virtual memory in use, in megabytes
      //! @return the virtual memory in use, in megabytes
      size_t GetVirtualMemoryInUse() const;

      //! @brief get the peak RAM used, in megabytes
      //! @return the peak RAM used, in megabytes
      size_t GetPeakRAMUsed() const;

      //! @brief get the RAM in use, in megabytes
      //! @return the RAM in use, in megabytes
      size_t GetRAMInUse() const;

      //! @brief write peak memory usage, if available, to an output stream
      //! @param STREAM stream to write the output to
      //! @return stream reference
      static std::ostream &WritePeakMemoryUsageInfo( std::ostream &STREAM);

      //! @brief write memory usage, if available, to an output stream
      //! @param STREAM stream to write the output to
      //! @return stream reference
      static std::ostream &WriteCurrentMemoryUsageInfo( std::ostream &STREAM);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MemoryUsage

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_MEMORY_USAGE_H_ 
