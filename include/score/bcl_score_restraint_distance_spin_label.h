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

#ifndef BCL_SCORE_RESTRAINT_DISTANCE_SPIN_LABEL_H_
#define BCL_SCORE_RESTRAINT_DISTANCE_SPIN_LABEL_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintDistanceSpinLabel
    //! @brief for scoring inter-CB distances in a protein model given a distance restraint.
    //!        It uses a histogram containing statistics for the frequency with which SL-CB distances are seen
    //!        over a database of proteins.
    //!
    //! @see @link example_score_restraint_distance_spin_label.cpp @endlink
    //! @author alexanns
    //! @date 12/14/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintDistanceSpinLabel :
      public RestraintAtomDistanceAssignment
    {
    private:

    //////////
    // data //
    //////////

      //! energy function with lower and upper bound
      static storage::Triplet< math::CubicSplineDamped, double, double> s_EnergyFunction;

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
      RestraintDistanceSpinLabel *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        return GetDefaultScheme();
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return GetDefaultScheme();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief operator() which takes an Assignment for calculating the agreement of the Assignment with the distance
      //! @param ASSIGNMENT the assignment which contains the MovableInterfaces whose distance will be scored
      //! @return return a double which is the score of the agreement of the MovableInterfaces with the distance
      double operator()
      (
        const restraint::AtomDistanceAssignment &ASSIGNMENT
      ) const;

      //! @brief ScoreDistance is the function which calculates the score a given distance should get
      //! @brief DISTANCE is the distance which will be scored - this is the distance from the protein model
      //! @brief RESTRAINT_DISTANCE is the distance that DISTANCE will be scored against
      //! @return returns the score of a distance of DISTANCE
      double ScoreDistance( const double &DISTANCE, const double &RESTRAINT_DISTANCE) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief gets the energy function and lower and upper bounds
      //! @return the energy function and lower and upper bounds
      static storage::Triplet< math::CubicSplineDamped, double, double> &GetEnergyFunction();

      //! @brief reads the energy function and lower and upper bounds
      //! @return the energy function and lower and upper bounds
      static storage::Triplet< math::CubicSplineDamped, double, double> ReadEnergyFunction();

    }; // class RestraintDistanceSpinLabel

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_DISTANCE_SPIN_LABEL_H_
