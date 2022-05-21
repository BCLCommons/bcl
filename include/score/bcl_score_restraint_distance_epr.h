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

#ifndef BCL_SCORE_RESTRAINT_DISTANCE_EPR_H_
#define BCL_SCORE_RESTRAINT_DISTANCE_EPR_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_bicubic_spline.h"
#include "restraint/bcl_restraint_atom_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintDistanceEPR
    //!
    //! @brief This class scores the agreement of a protein model with distance measurements obtained from an EPR
    //! experiment. The exposure of the spin labeling site is taken into account as well as the projection angles
    //! between the CaCb and CaCa vectors of the spin labeling site. By using these two parameters it is possible to
    //! determine if the spin labels are pointing towards or away from each other.
    //!
    //! @detail The neighbor count is used to quantify the exposure of a residue. For computing the neighbor count a
    //! minimum sequence separation of 0 and max euclidean distance of 11.4A are used. The exposure of the spin labeling
    //! site and the projection angles are aggregated with the formula arg = (nc_a * proj_a) + (nc_b * proj_b). The
    //! hereby computed value and the difference d between spin-spin and backbone distance are the arguments of the
    //! energy function E = f(arg, d).
    //!
    //! @see @link example_score_restraint_distance_epr.cpp @endlink
    //! @author fischea
    //! @date May 19, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintDistanceEPR :
      public ProteinModel
    {

    //////////
    // data //
    //////////

    private:

      //! scheme of this score
      std::string m_Scheme;

      //! shared pointer to the energy function which scores the agreement of a protein model with EPR distance restraints
      util::ShPtr< math::BicubicSpline> m_EnergyFunction;

      //! shared pointers to the EPR restraints
      util::ShPtrVector< restraint::AtomDistance> m_Restraints;

      //! shared pointer to the neighbor list calculator used to compute the exposure of the spin labeling sites
      util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > m_NeighborCalculator;

      //! neighbor count calculator used to quantify the exposure of the spinlabeling sites
      assemble::AANeighborCount m_ExposureCalculator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief construct from a restraint set, histogram file, and scheme
      //! @param RESTRAINTS EPR restraints for scoring the protein model
      //! @param HISTOGRAM_FILENAME name of the histogram file used to interpolate the potential function
      //! @param SCHEME scheme of this score
      RestraintDistanceEPR
      (
        const util::ShPtrVector< restraint::AtomDistance> &RESTRAINTS,
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief returns a pointer to a new RestraintDistanceEPR
      //! @return pointer to a new RestraintDistanceEPR
      RestraintDistanceEPR *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the scheme of this score
      //! @return the scheme of this score
      const std::string &GetScheme() const;

      //! @brief returns default file name where the statistics and in consequence the energy potentials are read from
      //! @return default file name where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns the default scheme of this score
      //! @return the default scheme of this score
      static const std::string &GetDefaultScheme();

    ////////////////
    // operations //
    ////////////////

      //! @brief scores the agreement of a given protein model with distance measurements obtained from an EPR experiment
      //! @param PROTEIN_MODEL protein model for which to compute the agreement
      //! @return normalized agreement score with -1 being the best and 0 being the worst agreement
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief reads in members from stream
      //! @param ISTREAM stream to read members from
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write members into a stream
      //! @param OSTREAM stream to write members into
      //! @INDENT number of indentations to use
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief creates the energy function used to score the agreement of model with the EPR data from the given histogram
      //! @param HISTOGRAM_FILENAME histogram from which to create the energy function
      //! @return shared pointer to the energy function created from the given histogram
      static util::ShPtr< math::BicubicSpline> CreateEnergyFunction( const std::string &HISTOGRAM_FILENAME);

    }; // class RestraintDistanceEPR

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_DISTANCE_EPR_H_
