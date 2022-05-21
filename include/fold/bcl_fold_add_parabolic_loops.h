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

#ifndef BCL_FOLD_ADD_PARABOLIC_LOOPS_H_
#define BCL_FOLD_ADD_PARABOLIC_LOOPS_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "linal/bcl_linal_vector_3d.h"
#include "restraint/bcl_restraint.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AddParabolicLoops
    //! @brief Add parabolic loop regions to protein models
    //!
    //! @see @link example_fold_add_parabolic_loops.cpp @endlink
    //! @author putnamdk, mendenjl
    //! @date Aug 28, 2012
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AddParabolicLoops :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! whether to determine the norm factor with regula falsi (true) or pythagorean approximation (false)
      bool m_DetermineAnalyticNormFactor;

      //! whether to approximate CB location
      bool m_ApproximateCB;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Constructor from members
      //! @param LOOPS bool value to represent loops that are not present in the protein model
      explicit AddParabolicLoops( const bool &ANALYTIC_NORM_FACTOR, const bool &APPROXIMATE_CB = false);

      //! @brief Clone function
      //! @return pointer to new AddParabolicLoops
      AddParabolicLoops *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief Get coordinates for residue placement of hypothetical loop
      //! @param PROTEIN_MODEL - protein model
      //! @param FIRST_RESIDUE, LAST_RESIDUE, first and last residues for the loop
      //! @return returns the coordinates for each atom in the loop
      storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> >
        GetLoopCoordinates( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief overloaded () operator to calculate Intensity from Q
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return returns Intensity for given Q for both the experimental and calculated data
      assemble::ProteinModel operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

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

    private:

      //! @brief return Computed Normalization factor for parabolic function - determines parabolic height.
      //! @param DESIRED_DISTANCE - this is how long you want the parabola to be
      //! @param SSE_DISTANCE - this is the head to tail distance between two sse's
      //! @return normalization factor
      double ComputeNormFactor( const double &DESIRED_DISTANCE, const double &SSE_DISTANCE) const;

      //! @brief return Computed Normalization factor for parabolic function - determines parabolic height.
      //! @param DESIRED_DISTANCE - this is how long you want the parabola to be
      //! @param SSE_DISTANCE - this is the head to tail distance between two sse's
      //! @return normalization factor
      double ComputeNormalizationFactor( const double &DESIRED_DISTANCE, const double &SSE_DISTANCE, const double &NAA) const;

    }; // class AddParabolicLoops

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_ADD_PARABOLIC_LOOPS_H_
