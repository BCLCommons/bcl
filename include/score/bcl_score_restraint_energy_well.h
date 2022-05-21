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

#ifndef BCL_SCORE_RESTRAINT_ENERGY_WELL_H_
#define BCL_SCORE_RESTRAINT_ENERGY_WELL_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_restraint_nmr_distance_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintEnergyWell
    //! @brief calculates a restraint score using an nergy well
    //! @details Calculates a Piecewise function from the
    //! experimental restraints already determined and determines a score based on the distances determined from the
    //! computed model Please see page 178 of the AMBER manual
    //! for R < r1, (with the slope of the "left-hand" parabola at the point R=r1) - 1
    //! for r1 <= R < r2, ( k2 ( R - r2 ) ^ 2) - 1
    //! for r2 <= R < r3, E = -1
    //! for r3 <= R < r4, ( k3 ( R - r3 ) ^ 2) -1
    //! for r4 <= R, ( with the slope of the "right-hand" parabola at the point R=r4) - 1
    //! where R is the model distance, r1 is some number, r2 is the lower experimental bound, r3 is the experimental
    //! distance, and r4 is the experimental upperbound
    //!
    //! @see @link example_score_restraint_energy_well.cpp @endlink
    //! @author akinlr, alexanns
    //! @date Jul 2, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintEnergyWell :
      public RestraintNMRDistanceInterface
    {

    public:

    //////////
    // enum //
    //////////

      //! enumerator for type of energy well
      enum Type
      {
        e_NOE,        //!< setup as described above
        e_PRE,        //!< no penalty when between lower and upper bounds, r3 is set to upper bound,
                      //!< with r1 and r4 given a set value away from the upper or lower bound
        s_NumberTypes //!< undefined type, also size of enum
      };

      //! @brief conversion to a string from a Type
      //! @param TYPE the type to get a string for
      //! @return a string representing that type
      static const std::string &GetTypeName( const Type &TYPE);

      //! @brief enum class wrapper for Type
      typedef util::WrapperEnum< Type, &GetTypeName, s_NumberTypes> TypeEnum;

    private:

      double m_KTwoValue;   //!< the K value to be used if the x value falls within the second range usually assigned 0
      double m_KThreeValue; //!< the K value to be used if the x value falls within the fourth range
      TypeEnum m_Type;          //!< the type of energy well to be used

      //! scheme to be used
      std::string m_Scheme;

      //! default r1 value
      static const double s_DefaultROne;

      //! distance between r1 and r2 as well as r3 and r4 for PRE calculations
      static const double s_PREPenaltyWidth;

      //! score for a restraint with residues/atoms not found in the protein model
      static const double s_DefaultScore;

      //! effective distance per bond
      static const double s_EffectiveDistancePerBond;

      //! instances for PRE and NOEs
      static const util::SiPtr< const util::ObjectInterface> s_PREInstance;
      static const util::SiPtr< const util::ObjectInterface> s_NOEInstance;

    public:

    //////////
    // data //
    //////////

      //! where the minima of the graph should fall
      static const double s_WellDepth;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintEnergyWell();

      //! @brief constructor from well type
      RestraintEnergyWell( const TypeEnum &WELL_TYPE);

      //! @brief parameter constructor
      //! @param K_TWO_VALUE the K value to be used if the x value falls within the second range
      //! @param K_THREE_VALUE the K value to be used if the x value falls within the fourth range
      //! @param WELL_TYPE type of energy well to be used
      //! @param SCHEME the short tag denoting this scoring function
      RestraintEnergyWell
      (
        const double K_TWO_VALUE,
        const double K_THREE_VALUE,
        const TypeEnum &WELL_TYPE = e_NOE,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief Clone function
      //! @return pointer to new RestraintEnergyWell
      virtual RestraintEnergyWell *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const;

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

      //! @brief create a PiecewiseFunction i.e. the scoring function specific to the current restraint distance
      //! @param RESTRAINT_DISTANCE gives the experimental data necessary for creating the function
      //! @param BOND_DISTANCE # of bonds from atoms in restraint to CB
      //! @return a PiecewiseFunction related to scoring NOEs
      math::PiecewiseFunction GetPiecewiseFunction
      (
        const restraint::AtomDistanceAssignment &RESTRAINT_DISTANCE,
        const size_t BOND_DISTANCE
      ) const;

      //! @brief will provide the Linear functions needed for the piecewise function because the slope is the same as
      //! the preceding quadratic function
      //! @param QUADRATIC gives the quadratic function needed to find the slope from its derivative
      //! @param R_VALUE gives the value needed to determine the y-intercept of the linear function
      //! @return a function interface in a linear function form
      static util::ShPtr< math::FunctionInterfaceSerializable< double, double> > GetLinearFunction
      (
        const math::QuadraticFunction &QUADRATIC, const double R_VALUE
      );

    }; // class RestraintEnergyWell

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_RESTRAINT_ENERGY_WELL_H_
