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

#ifndef BCL_SCORE_RESTRAINT_NOE_ATTRACTION_H_
#define BCL_SCORE_RESTRAINT_NOE_ATTRACTION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_restraint_atom_attraction.h"
#include "bcl_score_restraint_nmr_distance_interface.h"
#include "biol/bcl_biol_atom_types.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintNoeAttraction
    //! @brief used to add a "penalty" function to the NOE Knowledge based Potential scoring function based on the atom
    //!        distance given in the histogram, specifically for atom distances that are too large
    //!
    //! @see @link example_score_restraint_noe_attraction.cpp @endlink
    //! @author akinlr, weinerbe
    //! @date 06/30/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintNoeAttraction :
      public RestraintNMRDistanceInterface
    {
    private:

    //////////
    // data //
    //////////

      //! map of histogram data
      static storage::Map
      <
        storage::Pair< biol::AtomType, size_t>,           // type (CB, H, or HA) and number of bonds
        RestraintAtomAttraction
      > s_HistogramData;

      //! score for a restraint with residues/atoms not found in the protein model
      static const double s_DefaultScore;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone is the virtual copy constructor
      RestraintNoeAttraction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return GetDefaultScheme();
      }

      //! @brief get the histograms
      //! @return the histograms
      static const storage::Map
      <
        storage::Pair< biol::AtomType, size_t>,
        RestraintAtomAttraction
      > &GetNOEHistogram();

    ///////////////
    // operators //
    ///////////////

      //! @brief () operator scores protein model
      //! @param RESTRAINT restraint to be scored
      //! @return score
      double operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief reads in the histograms with NOE distance statistics from the file
      //! @return the histograms
      static storage::Map
      <
        storage::Pair< biol::AtomType, size_t>,
        RestraintAtomAttraction
      > ReadNOEHistogram();

    }; // class RestraintNoeAttraction
  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_NOE_ATTRACTION_H_
