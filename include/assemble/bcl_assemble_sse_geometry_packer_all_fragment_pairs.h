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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_ALL_FRAGMENT_PAIRS_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_ALL_FRAGMENT_PAIRS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_packer_interface.h"
#include "bcl_assemble_sse_geometry_packing.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackerAllFragmentPairs
    //! @brief class calculates the packing between all sub-geometry pairs for a given pair of geometries
    //! @details This class calculates all packings for all the sub-geometry pairs for a given pair of SSEGeometryInterface
    //! derived classes and returns the packings in a list of lists.
    //! For GEOMETRY_A ( composed of M sub-geometries) and GEOMETRY_B (composed of N sub-geometries)
    //! the returned list will have M list with each such list having N packings.
    //!
    //! @see @link example_assemble_sse_geometry_packer_all_fragment_pairs.cpp @endlink
    //! @author karakam
    //! @date Jul 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackerAllFragmentPairs :
      public SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > >
    {

    private:

    //////////
    // data //
    //////////

      //! minimal interface length
      double m_MinimalInterfaceLength;

      //! whether use distance cutoffs to decide which packings should be calculated
      bool m_UseDistanceCutoffs;

      //! map of distance cutoffs for each contact type
      storage::Map< biol::SSType, storage::Map< biol::SSType, double> > m_DistanceCutoffs;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEGeometryPackerAllFragmentPairs();

      //! @brief constructor from a minimal interface length
      //! @param MINIMAL_INTERFACE_LENGTH minimal interface length
      SSEGeometryPackerAllFragmentPairs
      (
        const double MINIMAL_INTERFACE_LENGTH
      );

      //! @brief constructor from a minimal interface length
      //! @param MINIMAL_INTERFACE_LENGTH minimal interface length
      //! @param DISTANCE_CUTOFFS map of distance cutoffs for each contact type
      SSEGeometryPackerAllFragmentPairs
      (
        const double MINIMAL_INTERFACE_LENGTH,
        const storage::Map< biol::SSType, storage::Map< biol::SSType, double> > &DISTANCE_CUTOFFS
      );

      //! @brief Clone function
      //! @return pointer to new SSEGeometryPackerAllFragmentPairs
      SSEGeometryPackerAllFragmentPairs *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the minimal interface length used
      //! @return minimal interface length used
      double GetMinimalInterfaceLength() const
      {
        return m_MinimalInterfaceLength;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the packing between all fragments of the given geometries and return it in a list of lists
      //! @param SSE_GEOMETRY_A first SSEGeometryInterface derived class of interest
      //! @param SSE_GEOMETRY_B second SSEGeometryInterface derived class of interest
      //! @return list of lists containing packing between all fragment pairings between given geometry pair
      storage::Vector< storage::List< SSEGeometryPacking> > operator()
      (
        const SSEGeometryInterface &SSE_GEOMETRY_A,
        const SSEGeometryInterface &SSE_GEOMETRY_B
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief return map of maps for upper distance ranges for sstype pairs
      //! @return map of maps for upper distance ranges for sstype pairs
      static const storage::Map< biol::SSType, storage::Map< biol::SSType, double> > &GetDistanceCutoffMap();

      //! @brief return map of maps for upper clash distance ranges for sstype pairs
      //! @return map of maps for upper clash distance ranges for sstype pairs
      static const storage::Map< biol::SSType, storage::Map< biol::SSType, double> > &GetClashDistanceCutoffMap();

    }; // class SSEGeometryPackerAllFragmentPairs

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKER_ALL_FRAGMENT_PAIRS_H_
