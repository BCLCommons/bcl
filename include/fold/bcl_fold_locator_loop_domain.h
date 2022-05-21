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

#ifndef BCL_FOLD_LOCATOR_LOOP_DOMAIN_H_
#define BCL_FOLD_LOCATOR_LOOP_DOMAIN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_locator_loop_segment.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "coord/bcl_coord_cyclic_coordinate_descent.h"
#include "find/bcl_find_locator_interface.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorLoopDomain
    //! @brief for locating and creating a loop domain in a protein model.
    //!
    //! @see @link example_fold_locator_loop_domain.cpp @endlink
    //! @author alexanns
    //! @date Aug 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorLoopDomain :
      public find::LocatorInterface< util::ShPtr< LoopDomain>, assemble::DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! the list of locators for the loop segments that make up the loop domain
      storage::List< LocatorLoopSegment> m_LoopSegments;

      //! true if the loop domain goes from n to c - false if LocatorLoopDomainCToN should be used instead
      bool m_NToCSequenceDirection;

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
      LocatorLoopDomain();

      //! @brief constructor taking member variable parameters
      //! @param SEGMENT_LOCATORS the loop segment locators that will be used to find the domain loop segments
      //! @param N_TO_C true if the loop domain goes from n to c - false if LocatorLoopDomainCToN should be used instead
      LocatorLoopDomain
      (
        const storage::List< LocatorLoopSegment> &SEGMENT_LOCATORS,
        const bool N_TO_C = true
      );

      //! @brief Clone function
      //! @return pointer to new LoopSegment
      LocatorLoopDomain *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief GetIdentification
      //! @return GetIdentification
      const std::string GetIdentification() const;

      //! @brief GetLoopSegments gives the list of locators for the loop segments that make up the loop domain
      //! @return the list of locators for the loop segments that make up the loop domain
      const storage::List< LocatorLoopSegment> &GetLoopSegments() const;

      //! @brief give vector of all the locator sses that correspond to this loop domain
      //! @return vector of all the locator sses that correspond to this loop domain
      util::SiPtrVector< const assemble::LocatorSSE> GetLocatorSSEs() const;

      //! @brief is this locator for the n to c seuqence direction
      //! @return true, if the locator is for c to n sequence direction (pseudo residue follows in sequence)
      bool IsNToCSequenceDirection() const
      {
        return m_NToCSequenceDirection;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate locates the desired LoopDomain in a domain
      //! @param SSE_DOMAIN the domain in which the LoopDomain will be located
      //! @return the desired loop domain as located in SSE_DOMAIN
      util::ShPtr< LoopDomain> Locate( const assemble::DomainInterface &SSE_DOMAIN) const;

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

      //! @brief calculates the sum of the square of the distances between target and moving points
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param LOCATOR_LOOP_DOMAIN Loop domain locator
      //! @param SUPERIMPOSE_ATOMS Set of atoms to be used in superimpostion
      //! @return calculates the sum of the square of the distances between target and moving points
      static double CalculateSquareDistanceSum
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const LocatorLoopDomain &LOCATOR_LOOP_DOMAIN
      );

      //! @brief calculates the sum of the square of the distances between target and moving points
      //! @param TARGET_AND_MOVING_POINTS List of target and moving points
      static double CalculateSquareDistanceSum
      (
        const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> &TARGET_AND_MOVING_POINTS
      );

      //! @brief calculates the root mean square distance between target and moving points
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param LOCATOR_LOOP_DOMAIN Loop domain locator
      //! @return the RMS of the points to be superimposed
      static double CalculateRMS
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const LocatorLoopDomain &LOCATOR_LOOP_DOMAIN
      );

      //! @brief determines if a loop domain is closed or not
      //! @param DOMAIN_LOCATOR the domain that will be checked for closure
      //! @param MODEL the protein model needed to see if loop is closed or not
      //! @param RMSD_CLOSURE_THRESHOLD the threshold in rmsd Angstrom for considering the loop closed
      //! @return bool true if the loop domain is closed - false otherwise
      static bool IsClosed
      (
        const LocatorLoopDomain &DOMAIN_LOCATOR,
        const assemble::ProteinModel &MODEL,
        const double RMSD_CLOSURE_THRESHOLD
      );

    }; // class LocatorLoopDomain

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOCATOR_LOOP_DOMAIN_H_ 
