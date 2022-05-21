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
#include "restraint/bcl_restraint_sas_optimization.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_saxs_opti_result.h"
#include "util/bcl_util_implementation.h"

using bcl::util::Message;

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

    const float SasOptimization::default_C1Min      = 0.80;
    const float SasOptimization::default_C1Max      = 1.20;
    const float SasOptimization::default_C2Min      = 0.00;
    const float SasOptimization::default_C2Max      = 4.00;
    const float SasOptimization::default_C1StepSize = 0.005;
    const float SasOptimization::default_C2StepSize = 0.100;

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasOptimization::s_Instance
    (
      GetObjectInstances().AddInstance( new SasOptimization())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasOptimization::SasOptimization() :
      m_c1Min( default_C1Min),
      m_c1Max( default_C1Max),
      m_c2Min( default_C2Min),
      m_c2Max( default_C2Max),
      m_c1StepSize( default_C1StepSize),
      m_c2StepSize( default_C2StepSize),
      m_scoreType( score::SasType::e_chi),
      m_useErrors( false),
      m_numEvaluations( 0),
      m_searchGrid( ( size_t) 41, ( size_t) 81, -100.0),
      m_ApproximateSideChains( false),
      m_Sans( false),
      m_DeuteriumExchangeParameter( 0.0)
    {

      storage::Vector< SasTransformation::TransformationTypeEnum> default_transform_vector;
      default_transform_vector.PushBack( SasTransformation::e_None);
      SasTransformation default_transform( default_transform_vector, false, m_useErrors, 1.0);

      m_transformation = default_transform;
    }

    //! @brief Constructor from members
    //! @param C1_MIN - minimum C1 value
    //! @param C1_MAX - maximum C1 value
    //! @param C2_MIN - minimum C2 value
    //! @param C2_MAX - maximum C2 value
    //! @param C1_STEPSIZE - step size for C1 parameter
    //! @param C2_STEPSIZE - step size for C2 parameter
    SasOptimization::SasOptimization
    (
      const float &C1_MIN,
      const float &C1_MAX,
      const float &C2_MIN,
      const float &C2_MAX,
      const float &C1_STEPSIZE,
      const float &C2_STEPSIZE,
      const score::SasType::ScoreFunctionEnum &SCORE_FUNCTION,
      const bool &USE_ERRORS,
      const SasTransformation &TRANSFORMATION_TYPES,
      const bool &APPROXIMATE_SIDE_CHAINS,
      const bool &HARDWARE_TYPE,
      const bool &SAS_TYPE,
      const float &DEUTERIUM_EXCHANGE_PARAMETER
    )
    :
      m_c1Min( C1_MIN),
      m_c1Max( C1_MAX),
      m_c2Min( C2_MIN),
      m_c2Max( C2_MAX),
      m_c1StepSize( C1_STEPSIZE),
      m_c2StepSize( C2_STEPSIZE),
      m_scoreType( SCORE_FUNCTION),
      m_useErrors( USE_ERRORS),
      m_numEvaluations( 0),
      m_searchGrid
      (
        (size_t)( (( m_c2Max - m_c2Min) / m_c2StepSize) + 1),
        (size_t)( (( m_c1Max - m_c1Min) / m_c1StepSize) + 1),
        ( float) - 100.0
      ),
      m_transformation( TRANSFORMATION_TYPES),
      m_ApproximateSideChains( APPROXIMATE_SIDE_CHAINS),
      m_Cpu( HARDWARE_TYPE),
      m_Sans( SAS_TYPE),
      m_DeuteriumExchangeParameter( DEUTERIUM_EXCHANGE_PARAMETER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasOptimization
    SasOptimization *SasOptimization::Clone() const
    {
      return new SasOptimization( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasOptimization::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief overloaded () operator to optimize c1 and c2 parameters
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      SaxsOptiResult result( GridWalk( PROTEIN_MODEL));

      BCL_MessageStd( "Transformation values:" + util::Format()( m_transformation));

      return result;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns lowest chi score from adusting c1 and c2 parameters over a grid of 3200 points
    //! @brief searches the grid by a line search using quadratic interpolation along each row of the grid
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::GridWalk( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      storage::Vector< float> chi_array;
      storage::Vector< float> c1_array;
      storage::Vector< float> c2_array;

      m_numEvaluations = 0;

      for( float c2( m_c2Min); c2 <= m_c2Max; c2 += m_c2StepSize)
      {
        SaxsOptiResult result( GoldenSection( m_c1Min, m_c1Max, c2, PROTEIN_MODEL));
        //SaxsOptiResult result( QuadraticInterpolation( m_c1Min, m_c1Max, c2, PROTEIN_MODEL));
        c1_array.PushBack( result.GetC1());
        c2_array.PushBack( result.GetC2());
        chi_array.PushBack( result.GetChi());
      }

      BCL_MessageStd( "chi_array: " + util::Format()( chi_array));

      //Get index of smallest element in list
      size_t index( GetMinimumIndex( chi_array));

      float final_chi( chi_array( index));
      float final_c1( c1_array( index));
      float final_c2( c2_array( index));

      SaxsOptiResult result( final_c1, final_c2, final_chi);

      return result;
    }

    //! @brief returns lowest chi score from adusting c1 and c2 parameters over a grid of 3200 points
    //! @brief for each column of the grid each row is tested until the function at a given point increases
    //! @brief this method is faster than a complete grid search, but still slow to find a final answer
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::SmartGridSearch( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize variables for grid search
      float chi( 0.0);
      float lowest_chi( 10000.0);
      float lowest_c1( 0.0);
      float lowest_c2( 0.0);

      float previous( 10000.0);
      float current( 0.0);
      size_t count( 0);

      // instead of running 1 calculation, I want to run a grid search on the c1 and c2 parameters
      // c1 = [0.80 - 1.20] step = 0.005   80 total points
      // c2 = [0.0 -  4.0] step = 0.1      40 total points
      // 3200 or less total evaluations for the optimimal c1 and c2 combination to minimize chi

      for( float c1( m_c1Min); c1 <= m_c1Max; c1 += m_c1StepSize)
      {
        for( float c2( m_c2Min); c2 <= m_c2Max; c2 += m_c2StepSize)
        {
          chi = EvaluateFunction( c1, c2, PROTEIN_MODEL);

          if( chi < lowest_chi)
          {
            lowest_chi = chi;
            lowest_c1 = c1;
            lowest_c2 = c2;
          }

          current = chi;
          if( current < previous)
          {
            previous = current;
          }
          else
          {
            previous = 10000.0;
            current = 0.0;
            break;
          }
          //BCL_Message( util::Message::e_Standard, "iteration: " + util::Format()( count));
          ++count;
        }
      }

      SaxsOptiResult result( lowest_c1, lowest_c2, lowest_chi);
      return result;
    }

    //! @brief one Dimensional unimodal derivative free optimization method using Quadratic Interpolation
    //! @param MIN_C1 - Minimum C1 value
    //! @param MAX_C1 - Maximum C1 value
    //! @param C2_VALUE - Constant row of the matrix
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::QuadraticInterpolation
    (
      float MIN_C1,
      float MAX_C1,
      float C2_VALUE,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      float a, b, c, d;
      float function_at_a, function_at_b, function_at_c, function_at_d;

      a = MIN_C1;
      b = MAX_C1;
      c = (MIN_C1 + MAX_C1) / 2;

      function_at_a = EvaluateFunction( a, C2_VALUE, PROTEIN_MODEL);
      function_at_b = EvaluateFunction( b, C2_VALUE, PROTEIN_MODEL);
      function_at_c = EvaluateFunction( c, C2_VALUE, PROTEIN_MODEL);

      d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

      // BCL_Message( util::Message::e_Standard, " d is: " + util::Format()( d));

      // handle the case where d is outside of the interval [a, b]
      if( !InsideQuadraticInterval( a, b, d))
      {
        SaxsOptiResult result;

        if( function_at_a < function_at_b)
        {
          result.SetC1( a);
          result.SetC2( C2_VALUE);
          result.SetChi( function_at_a);
        }
        else
        {
          result.SetC1( b);
          result.SetC2( C2_VALUE);
          result.SetChi( function_at_b);
        }
        return result;
      }

      // handle the case where c = d
      if( c == d)
      {
        c = c + 0.001;
      }

      function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);

      float prev_d( d);
      do
      {
        prev_d = d;
        // f increases on [C,D]
        if( c < d && function_at_c < function_at_d)
        {
          // a = a
          // c = c
          b = d;
          function_at_b = function_at_d;
          d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);
        }

        // f decreases on [D,C]
        else if( c > d && function_at_c < function_at_d)
        {
          a = d;
          function_at_a = function_at_d;
          // c = c
          // b = b;
          d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);
        }
        else if( c < d && function_at_c > function_at_d)
        {
          a = c;
          function_at_a = function_at_c;

          c = d;
          function_at_c = function_at_d;

          // b = b;
          d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);
        }
        else
        {
          // a = a;

          b = c;
          function_at_b = function_at_c;

          c = d;
          function_at_c = function_at_d;

          // b = b;
          d = Quadratic( a, function_at_a, b, function_at_b, c, function_at_c);

          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);
        }

      } while( InsideQuadraticInterval( MIN_C1, MAX_C1, d) && ( d - prev_d >= 0.002));

      SaxsOptiResult result;

      if( function_at_c < function_at_d)
      {
        result.SetC1( c);
        result.SetC2( C2_VALUE);
        result.SetChi( function_at_c);
      }
      else
      {
        result.SetC1( d);
        result.SetC2( C2_VALUE);
        result.SetChi( function_at_d);

      }
      return result;
    }

    //! @brief one Dimensional unimodal derivative free optimization method using GoldenSection Method
    //! @param MIN_C1 - Minimum C1 value
    //! @param MAX_C1 - Maximum C1 value
    //! @param C2_VALUE - Constant row of the matrix
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return returns storage vector of 3 values: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::GoldenSection
    (
      float MIN_C1,
      float MAX_C1,
      float C2_VALUE,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // This is an implementation of the Golden Section algorithm
      // first initialize all variables;

      // golden ratio (sqrt(5) - 1) / 2
      float theta( 0.618033989);

      float a( MIN_C1);
      float function_at_a( EvaluateFunction( a, C2_VALUE, PROTEIN_MODEL));

      float b( MAX_C1);
      float function_at_b( EvaluateFunction( b, C2_VALUE, PROTEIN_MODEL));

      float c( theta * a + ( 1 - theta) * b);
      float function_at_c( EvaluateFunction( c, C2_VALUE, PROTEIN_MODEL));

      float d( ( 1 - theta) * a + theta * b);
      float function_at_d( EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL));

      storage::Vector< float> storage_function_values
      (
        storage::Vector< float>::Create( function_at_a, function_at_b, function_at_c, function_at_d)
      );

      float prev_min_value( GetMinimumValue( storage_function_values));
      float current_min_value( prev_min_value);

      do
      {
        prev_min_value = current_min_value;

        if( function_at_c < function_at_d)
        {
          // a remains the same
          b = d;
          function_at_b = function_at_d;

          d = c;
          function_at_d = function_at_c;

          c =    theta * a + ( 1 - theta)  * b;
          function_at_c = EvaluateFunction( c, C2_VALUE, PROTEIN_MODEL);

          storage::Vector< float> c_less_than_d
          (
            storage::Vector< float>::Create( function_at_a, function_at_b, function_at_c, function_at_d)
          );

          current_min_value = GetMinimumValue( c_less_than_d);
        }
        else
        {
          // b remains the same
          a = c;
          function_at_a = function_at_c;

          c = d;
          function_at_c = function_at_d;

          d = ( 1 - theta) * a + theta * b;
          function_at_d = EvaluateFunction( d, C2_VALUE, PROTEIN_MODEL);

          storage::Vector< float> c_greater_equal_d
          (
            storage::Vector< float>::Create( function_at_a, function_at_b, function_at_c, function_at_d)
          );

          current_min_value = GetMinimumValue( c_greater_equal_d);
        }
      } while( math::Absolute( double( c - a)) > ( 2 * m_c1StepSize));

      //while ( math::Absolute( double( current_min_value - prev_min_value)) > m_c1StepSize);

      storage::Vector< SaxsOptiResult> result_vector;
      storage::Vector< float> values;

      SaxsOptiResult result_a( a, C2_VALUE, function_at_a);
      result_vector.PushBack( result_a);
      values.PushBack( function_at_a);

      SaxsOptiResult result_b( b, C2_VALUE, function_at_b);
      result_vector.PushBack( result_b);
      values.PushBack( function_at_b);

      SaxsOptiResult result_c( c, C2_VALUE, function_at_c);
      result_vector.PushBack( result_c);
      values.PushBack( function_at_c);

      SaxsOptiResult result_d( d, C2_VALUE, function_at_d);
      result_vector.PushBack( result_d);
      values.PushBack( function_at_d);

      size_t index( GetMinimumIndex( values));

      return result_vector( index);
    }

    //! @brief returns lowest chi score from a search pattern starting at quadratic interpolation line search
    //! @param RADIUS, search radius
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return SaxsOptiResult: c1, c2, and min chi score
    SaxsOptiResult SasOptimization::HookJeeves
    (
      linal::VectorND< float, 2> INITIAL_POINT,
      float RADIUS,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      //restraint::SaxsOptiResult start_point( QuadraticInterpolation( m_c1Min, m_c1Max, m_c2Min, PROTEIN_MODEL));

      linal::VectorND< float, 2> grid_point( DiscretePoint( INITIAL_POINT));

      storage::Vector< linal::VectorND< float, 2> > x_coordinates;
      x_coordinates.PushBack( grid_point);
      linal::VectorND< float, 2> x;
      float h( RADIUS);
      size_t k( 0);

      float epsilon( 1.0);
      while( h >= epsilon)
      {
        bool found( false);

        while( found == false && h >= epsilon)
        {
          x = Explore( x_coordinates( k), h, PROTEIN_MODEL);
          if( x != x_coordinates( k))
          {
            x_coordinates.PushBack( x);
            k = k + 1;
            found = true;
          }
          else
          {
            h = h / 2;
          }
        }
      }

      linal::VectorND< float, 2> final_x = x_coordinates.LastElement();

      float min_value( EvaluatePoint( final_x, PROTEIN_MODEL));

      SaxsOptiResult result( final_x.First(), final_x.Last(), min_value);
      return result;
    }

    // helper functions

    //! @brief Convert C1 coordinates to Matrix position
    //! @param C1 - value of contrast density parameter
    //! @return matrix location of C1 value
    int SasOptimization::C1ParamtoMatrix( float C1) const
    {
      float x( ( C1 - m_c1Min) / m_c1StepSize);
      return x > 0 ? ( x + 0.5) : ( x - 0.5);
    }

    //! @brief Convert C2 coordinates to Matrix position
    //! @param C2 - value of contrast density parameter
    //! @return matrix location of C2 value
    int SasOptimization::C2ParamtoMatrix( float C2) const
    {
      float x( ( C2 - m_c2Min) / m_c2StepSize);
      return x > 0 ? ( x + 0.5) : ( x - 0.5);
    }

    //! @brief takes input point from line search method and rounds it to a valid start point on the grid
    //! @param INITIAL_POINT point from line search that is not on a grid point
    //! @return point on Grid
    linal::VectorND< float, 2> SasOptimization::DiscretePoint( const linal::VectorND< float, 2> INITIAL_POINT)const
    {
      float c1( INITIAL_POINT.First());
      float c2( INITIAL_POINT.Last());

      int c1_int( c1 * 100);
      int c2_int( c2 * 10);

      float c1_grid( (float)c1_int / 100);
      float c2_grid( (float)c2_int / 10);
      linal::VectorND< float, 2> final_point;
      final_point(0) = c1_grid;
      final_point(1) = c2_grid;
      return final_point;
    }

    //! @brief compute chi for given parameters and protein model
    //! @param C1 - electron density adjustment parameter for form factors
    //! @param C2 - hydration shell adjustmet parameter for form factors
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return chi score between experimental and computed saxs profiles
    float SasOptimization::EvaluateFunction( float C1, float C2, const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      m_numEvaluations++;
      BCL_Message
      (
        util::Message::e_Standard,
           " evaluation number: " + util::Format()( m_numEvaluations)
      );

      // Setup Commandline Strings for either the opencl or non-opencl version of the code
      util::Implementation< SasDebyeInterface> sas
      (
        SasAnalysis::SetDebyeImplementation
        ( false, false, C1, C2, m_Cpu, m_Sans, m_DeuteriumExchangeParameter)
      );

      sas->SetExperimentalData( m_ExpData);

      SasExperimentalAndCalculatedData data_sets( sas->operator()( PROTEIN_MODEL));
      SasExperimentalAndCalculatedData transformed_data( m_transformation( data_sets));

      float score = ( float) score::SasType( m_useErrors, m_scoreType)( transformed_data);

      BCL_MessageStd( "C1: " + util::Format()( C1) + " C2: " + util::Format()( C2) + " Score: " + util::Format()( score));

      return score;
    }

    //! @brief returns lowest chi score all points around the center point based on Radius
    //! @param POINT, center point consisting of c1 and c2
    //! @param RADIUS, search radius
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return linal::VectorND< float, 2>: c1, c2 coordinates of min function value
    linal::VectorND< float, 2> SasOptimization::Explore
    (
      const linal::VectorND< float, 2> &POINT,
      float RADIUS,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // These are the direction unit vectors
      linal::VectorND< float, 2> e1;
      e1( 0) = m_c1StepSize;
      e1( 1) = 0;

      linal::VectorND< float, 2> e2;
      e2( 0) = 0;
      e2( 1) = m_c2StepSize;

      // These are the diagonal directions
      linal::VectorND< float, 2> e3;
      e3( 0) = m_c1StepSize;
      e3( 1) = m_c2StepSize;

      linal::VectorND< float, 2> e4;
      e4( 0) = -m_c1StepSize;
      e4( 1) = -m_c2StepSize;

      linal::VectorND< float, 2> e5;
      e5( 0) = m_c1StepSize;
      e5( 1) = -m_c2StepSize;

      linal::VectorND< float, 2> e6;
      e6( 0) = -m_c1StepSize;
      e6( 1) = m_c2StepSize;

      float h( RADIUS);

      linal::VectorND< float, 2> x0 = POINT;

      // pick initial evaluation points
      linal::VectorND< float, 2> x2, x3, x4, x5, x6, x7, x8, x9;
      x2 = x0 + h * e1;
      x3 = x0 + h * e2;
      x4 = x0 - h * e1;
      x5 = x0 - h * e2;

      x6 = x0 + h * e3;
      x7 = x0 + h * e4;
      x8 = x0 + h * e5;
      x9 = x0 + h * e6;

      storage::Vector< linal::VectorND< float, 2> > search_points;
      search_points.PushBack( x0);
      search_points.PushBack( x2);
      search_points.PushBack( x3);
      search_points.PushBack( x4);
      search_points.PushBack( x5);
      search_points.PushBack( x6);
      search_points.PushBack( x7);
      search_points.PushBack( x8);
      search_points.PushBack( x9);

      storage::Vector< float> function_values;

      for
      (
        storage::Vector< linal::VectorND< float, 2> >::const_iterator
          itr_search_points( search_points.Begin()),
          itr_search_points_end( search_points.End());
        itr_search_points != itr_search_points_end;
        ++itr_search_points
      )
      {
        if( Inbounds( *itr_search_points))
        {
          float value = EvaluatePoint( *itr_search_points, PROTEIN_MODEL);
          function_values.PushBack( value);
        }
        else
        {
          function_values.PushBack( 99999);
        }
      }

      size_t min_index( GetMinimumIndex( function_values));

      return search_points( min_index);
    }

    //! @brief returns the function value of a given c1 and c2 combination.  Results are stored in m_searchGrid to
    //! @brief prevent duplication of work
    //! @param POINT, center point consisting of c1 and c2
    //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
    //! @return function value of given c1 and c2 parameter
    float SasOptimization::EvaluatePoint
    (
      const linal::VectorND< float, 2> &POINT,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      float function_at_x = util::GetUndefined< float>();

      float c1( POINT.First());
      float c2( POINT.Last());

      int c1_matrix = C1ParamtoMatrix( c1);
      int c2_matrix = C2ParamtoMatrix( c2);

      if( m_searchGrid( c2_matrix, c1_matrix) == -100.0)
      {
        function_at_x = EvaluateFunction( c1, c2, PROTEIN_MODEL);
        m_searchGrid( c2_matrix, c1_matrix) = function_at_x;
      }
      else
      {
        function_at_x = m_searchGrid( c2_matrix, c1_matrix);
      }

      return function_at_x;
    }

    //! @brief find index of minimum value in storage vector
    //! @param ARRAY storage vector to search
    //! @return index of minimum value in array
    size_t SasOptimization::GetMinimumIndex( const storage::Vector< float> &ARRAY) const
    {
      size_t size( ARRAY.GetSize());
      size_t index( 0);

      for( size_t i( 1); i < size; ++i)
      {
        if( ARRAY( i) < ARRAY( index))
        {
          index = i;
        }
      }

      return index;
    }

    //! @brief return minimum value in storage vector
    //! @param ARRAY storage vector to search
    //! @return minimum value in array
    float SasOptimization::GetMinimumValue( const storage::Vector< float> &ARRAY) const
    {
      linal::Vector< float> linal_function_values( ARRAY);
      return linal_function_values.Min();
    }

    //! @brief Ensure that the proposed point to evaluate is inside the searchGrid
    //! @param POINT, center point consisting of c1 and c2
    //! @return true if inbounds
    bool SasOptimization::Inbounds( linal::VectorND< float, 2> POINT) const
    {
      bool inbounds( true);

      if( POINT.First() < m_c1Min)
      {
        inbounds = false;
      }
      else if( POINT.First() > m_c1Max)
      {
        inbounds = false;
      }

      else if( POINT.Last() < m_c2Min)
      {
        inbounds = false;
      }
      else if( POINT.Last() > m_c2Max)
      {
        inbounds = false;
      }

      return inbounds;
    }

    //! @brief Ensure that the proposed point to evaluate (D) is inside the specified interval
    //! @param A - Left Boundary Point
    //! @param B - Right Boundary Point
    //! @param TEST_POINT, point to evaluate
    //! @return true if inbounds
    bool SasOptimization::InsideQuadraticInterval( float A, float B, float TEST_POINT) const
    {
      bool inbounds( true);
      if( TEST_POINT < A)
      {
        inbounds = false;
      }
      else if( TEST_POINT > B)
      {
        inbounds = false;
      }
      return inbounds;
    }

    //! @brief computes a quadratic model based on 3 values and thier function
    //! @param A - First point
    //! @param FA - function at first point
    //! @param B - Second point
    //! @param FB - function at second point
    //! @param C - Third point
    //! @param FC - function at third point
    float SasOptimization::Quadratic( float A, float FA, float B, float FB, float C, float FC) const
    {
      float calculation
      (
        0.5
        * ( FA * ( B*B - C*C ) - FC * ( B*B - A*A) + FB * ( C*C - A*A))
        / ( FA * ( B - C)      - FC * ( B - A)     + FB * ( C - A))
       );

      return calculation;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasOptimization::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_c1Min, ISTREAM);
      io::Serialize::Read( m_c1Max, ISTREAM);
      io::Serialize::Read( m_c2Min, ISTREAM);
      io::Serialize::Read( m_c2Max, ISTREAM);
      io::Serialize::Read( m_c1StepSize, ISTREAM);
      io::Serialize::Read( m_c2StepSize, ISTREAM);
      io::Serialize::Read( m_numEvaluations, ISTREAM);
      io::Serialize::Read( m_searchGrid, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasOptimization::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_c1Min, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c1Max, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c2Min, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c2Max, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c1StepSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_c2StepSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_numEvaluations, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_searchGrid, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
