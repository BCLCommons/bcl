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
#include "score/bcl_score_restraint_energy_well.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_const_function.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_quadratic_function.h"
#include "math/bcl_math_trigonometric_transition.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // enum //
  //////////

    //! @brief conversion to a string from a Type
    //! @param TYPE the type to get a string for
    //! @return a string representing that type
    const std::string &RestraintEnergyWell::GetTypeName( const RestraintEnergyWell::Type &TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "NOE",
        "PRE",
        GetStaticClassName< Type>()
      };
      return s_descriptors[ size_t( TYPE)];
    }

  //////////
  // data //
  //////////

    //! default r1 value
    const double RestraintEnergyWell::s_DefaultROne( 0.0);

    //! distance between r1 and r2 as well as r3 and r4 for PRE calculations
    const double RestraintEnergyWell::s_PREPenaltyWidth( 15.0);

    //! initialize "s_WellDepth"
    const double RestraintEnergyWell::s_WellDepth( -1.0);

    //! score for a restraint with residues/atoms not found in the protein model
    const double RestraintEnergyWell::s_DefaultScore( 0.0);

    //! effective distance per bond
    const double RestraintEnergyWell::s_EffectiveDistancePerBond( 1.0);

    const util::SiPtr< const util::ObjectInterface> RestraintEnergyWell::s_PREInstance
    (
      util::Enumerated< RestraintAtomDistanceAssignment>::AddInstance( new RestraintEnergyWell( e_PRE))
    );
    const util::SiPtr< const util::ObjectInterface> RestraintEnergyWell::s_NOEInstance
    (
      util::Enumerated< RestraintAtomDistanceAssignment>::AddInstance( new RestraintEnergyWell( e_NOE))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintEnergyWell::RestraintEnergyWell() :
      m_KTwoValue( util::GetUndefinedDouble()),
      m_KThreeValue( util::GetUndefinedDouble()),
      m_Type( s_NumberTypes),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief constructor from well type
    RestraintEnergyWell::RestraintEnergyWell( const TypeEnum &WELL_TYPE) :
      m_KTwoValue( 1.0),
      m_KThreeValue( 1.0),
      m_Type( WELL_TYPE),
      m_Scheme( GetDefaultScheme() + "_" + WELL_TYPE.GetString())
    {
    }

    //! @brief parameter constructor
    //! @param K_TWO_VALUE the K value to be used if the x value falls within the second range
    //! @param K_THREE_VALUE the K value to be used if the x value falls within the fourth range
    //! @param WELL_TYPE type of energy well to be used
    //! @param SCHEME the short tag denoting this scoring function
    RestraintEnergyWell::RestraintEnergyWell
    (
      const double K_TWO_VALUE,
      const double K_THREE_VALUE,
      const TypeEnum &WELL_TYPE,
      const std::string &SCHEME
    ) :
      m_KTwoValue( K_TWO_VALUE),
      m_KThreeValue( K_THREE_VALUE),
      m_Type( WELL_TYPE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RestraintEnergyWell
    RestraintEnergyWell *RestraintEnergyWell::Clone() const
    {
      return new RestraintEnergyWell( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RestraintEnergyWell::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &RestraintEnergyWell::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "energy_well");

      // end
      return s_default_scheme;
    }

    //! @brief returns scheme being used
    //! @return scheme being used
    const std::string &RestraintEnergyWell::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief () operator scores protein model
    //! @param RESTRAINT restraint to be scored
    //! @return score
    double RestraintEnergyWell::operator()( const restraint::AtomDistanceAssignment &RESTRAINT) const
    {
      // calculate the distance
      const double cb_distance( RESTRAINT.CalculateAtomDistance());

      // if the calculated distance is undefined
      if( !util::IsDefined( cb_distance))
      {
        // return default score
        return s_DefaultScore;
      }

      // get the bond distance from the CB
      const size_t bond_distance( GetTotalBondsFromCB( RESTRAINT));

      // if the bond distance is not defined
      if( !util::IsDefined( bond_distance))
      {
        // return default score
        return s_DefaultScore;
      }

      // create a Piecewise Function
      const math::PiecewiseFunction noe_score_piecewise_function
      (
        GetPiecewiseFunction( RESTRAINT, bond_distance)
      );

      // return the piecewise function with the calculated NOE distance restraint score from BCL data -1 in order to
      // align it with the other three scoring functions
      return noe_score_piecewise_function( cb_distance) + s_WellDepth;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create a PiecewiseFunction i.e. the scoring function specific to the current restraint distance
    //! @param RESTRAINT_DISTANCE gives the experimental data necessary for creating the function
    //! @param BOND_DISTANCE # of bonds from atoms in restraint to CB
    //! @return a PiecewiseFunction related to scoring NOEs
    math::PiecewiseFunction RestraintEnergyWell::GetPiecewiseFunction
    (
      const restraint::AtomDistanceAssignment &RESTRAINT_DISTANCE,
      const size_t BOND_DISTANCE
    ) const
    {
      // make sure the restraint type is defined
      BCL_Assert( m_Type != s_NumberTypes, "Energy well type must be defined");

      // initialize r-values
      double r_one;
      double r_two;
      double r_three;
      double r_four;

      // calculate the effective distance
      const double effective_distance( double( BOND_DISTANCE) * s_EffectiveDistancePerBond);

      // if the type is NOE
      if( m_Type == e_NOE)
      {
        // pull each necessary element from the parameters of ScoreDistance
        const double upper_bound( RESTRAINT_DISTANCE.GetUpperBound() + effective_distance);
        const double lower_bound( std::max( s_DefaultROne, RESTRAINT_DISTANCE.GetLowerBound() - effective_distance));
        const double exp_distance( RESTRAINT_DISTANCE.GetDistance() + effective_distance);

        r_one = std::min( s_DefaultROne, lower_bound - 0.1); // this is to make sure r_one and r_two are not the same
        r_two = lower_bound;
        r_three = exp_distance;
        r_four = std::max( exp_distance + 0.1, upper_bound); // this is to make sure r_three and r_four are not the same
      }
      // if the type is PRE
      else // if( m_Type == e_PRE)
      {
        // shift the bounds since the SL distance is usually longer than the CB distance
        const double upper_bound( RESTRAINT_DISTANCE.GetUpperBound() - effective_distance / 2.0);
        const double lower_bound
        (
          std::max( s_DefaultROne, RESTRAINT_DISTANCE.GetLowerBound() - effective_distance / 2.0)
        );

        r_one = lower_bound - s_PREPenaltyWidth;
        r_two = lower_bound;
        r_three = upper_bound;
        r_four = upper_bound + s_PREPenaltyWidth;
      }

      const double max( std::numeric_limits< double>::max());
      // assert that the r-values are consecutive
      BCL_Assert
      (
        r_one < r_two && r_two < r_three && r_three < r_four,
        "Energy well score r-values must be sequential: " + util::Format()( r_one) + ", " +
        util::Format()( r_two) + ", " + util::Format()( r_three) + ", " + util::Format()( r_four)
      );

      // use these parameters to construct the ranges necessary for the piecewise function
      const math::Range< double>
        range_a( math::RangeBorders::e_LeftOpen, -max, r_one, math::RangeBorders::e_RightOpen),
        range_b( math::RangeBorders::e_LeftClosed, r_one, r_two, math::RangeBorders::e_RightOpen),
        range_c( math::RangeBorders::e_LeftClosed, r_two, r_three, math::RangeBorders::e_RightOpen),
        range_d( math::RangeBorders::e_LeftClosed, r_three, r_four, math::RangeBorders::e_RightOpen),
        range_e( math::RangeBorders::e_LeftClosed, r_four, max, math::RangeBorders::e_RightOpen);

      // create the vectors needed for the Quadratic Functions
      // these are the x and y values of vertex of the two quadratic functions
      const storage::VectorND< 2, double> bound_two( r_two, 0.0);
      const storage::VectorND< 2, double> bound_three( r_three, 0.0);

      // create the functions
      const util::ShPtr< math::FunctionInterfaceSerializable< double, double> > function_c
      (
        new math::ConstFunction< double, double>( 0.0)
      );
      util::ShPtr< math::FunctionInterfaceSerializable< double, double> > function_a;
      util::ShPtr< math::FunctionInterfaceSerializable< double, double> > function_b;
      util::ShPtr< math::FunctionInterfaceSerializable< double, double> > function_d;
      util::ShPtr< math::FunctionInterfaceSerializable< double, double> > function_e;

      // if NOE
      if( m_Type == e_NOE)
      {
        // create Quadratic Functions so the Linear functions can be calculated
        const math::QuadraticFunction quad_a( bound_two, m_KTwoValue);
        const math::QuadraticFunction quad_e( bound_three, m_KThreeValue);

        // use parameters from ScoreDistance to construct the functions necessary for the piecewise function
        function_a = GetLinearFunction( quad_a, r_one);
        function_b = util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
          (
            new math::QuadraticFunction( bound_two, m_KTwoValue)
          );
        function_d = util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
          (
            new math::QuadraticFunction( bound_three, m_KThreeValue)
          );
        function_e = GetLinearFunction( quad_e, r_four);
      }
      // PRE
      else
      {
        // use cosine transition instead of quadratic
        function_a = util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
          (
            new math::ConstFunction< double, double>( 1.0)
          );
        function_b = util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
          (
            new math::TrigonometricTransition( r_one, bound_two.First(), 1.0, bound_two.Second())
          );
        function_d = util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
          (
            new math::TrigonometricTransition( bound_three.First(), r_four, bound_three.Second(), 1.0)
          );
        function_e = function_a;
      }

      // put the ranges and corresponding functions into pairs into a list
      storage::List
      <
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >
      > noe_score_list;
      noe_score_list.PushBack
      (
        storage::Pair
        <
          math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
        >( range_a, function_a)
      ),
      noe_score_list.PushBack
      (
        storage::Pair
        <
          math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
        >( range_b, function_b)
      ),
      noe_score_list.PushBack
      (
        storage::Pair
        <
          math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
        >( range_c, function_c)
      ),
      noe_score_list.PushBack
      (
        storage::Pair
        <
          math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
        >( range_d, function_d)
      ),
      noe_score_list.PushBack
      (
        storage::Pair
        <
          math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> >
        >( range_e, function_e)
      );

      // construct the piecewise function
      return math::PiecewiseFunction( noe_score_list);
    }

    //! @brief will provide the Linear functions needed for the piecewise function because the slope is the same as the
    //! preceding quadratic function
    //! @param QUADRATIC gives the quadratic function needed to find the slope from its derivative
    //! @param R_VALUE gives the value needed to determine the y-intercept of the linear function
    //! @return a function interface in a linear function form
    util::ShPtr< math::FunctionInterfaceSerializable< double, double> > RestraintEnergyWell::GetLinearFunction
    (
      const math::QuadraticFunction &QUADRATIC, const double R_VALUE
    )
    {
      // get the derivative of "QUADRATIC"
      const util::ShPtr< math::FunctionInterfaceSerializable< double, double> > derivative( QUADRATIC.GetDerivative());

      // use "derivative" to get the slope of "QUADRATIC" at "R_VALUE"
      const double slope( derivative->operator()( R_VALUE));

      // get the y_value of "QUADRATIC" at "R_VALUE"
      const double y_value( QUADRATIC( R_VALUE));

      // solve for the y-intercept and return the corresponding LinearFunction
      const double y_intercept( y_value - ( slope * R_VALUE));
      return util::ShPtr< math::FunctionInterfaceSerializable< double, double> >( new math::LinearFunction( slope, y_intercept));
    }

  } // namespace score
} // namespace bcl
