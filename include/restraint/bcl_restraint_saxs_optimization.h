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

#ifndef BCL_RESTRAINT_SAXS_OPTIMIZATION_H_
#define BCL_RESTRAINT_SAXS_OPTIMIZATION_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_saxs_transformation.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector_nd.h"
#include "score/bcl_score_saxs_type.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SaxsOptimization
    //! @brief Class to optimize c1 and c2 parameters
    //!
    //! @see @link example_restraint_SaxsOptimization.cpp @endlink
    //! @author putnamdk
    //! @date May 29, 2013
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SaxsOptimization :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      float m_c1Min;
      float m_c1Max;
      float m_c2Min;
      float m_c2Max;
      float m_c1StepSize;
      float m_c2StepSize;
      score::SaxsType::ScoreFunctionEnum m_scoreType;
      bool m_useErrors;
      mutable size_t m_numEvaluations;
      mutable linal::Matrix< float> m_searchGrid;
      SaxsTransformation m_transformation;
      bool m_ApproximateSideChains;
      bool m_Cpu;

      //! Experimental data
      util::ShPtr< SaxsScatteringData> m_ExpData;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SaxsOptimization();

      //! @brief Constructor from members
      //! @param C1_MIN - minimum C1 value
      //! @param C1_MAX - maximum C1 value
      //! @param C2_MIN - minimum C2 value
      //! @param C2_MAX - maximum C2 value
      //! @param C1_STEPSIZE - step size for C1 parameter
      //! @param C2_STEPSIZE - step size for C2 parameter
      SaxsOptimization
      (
        const float &C1_MIN,
        const float &C1_MAX,
        const float &C2_MIN,
        const float &C2_MAX,
        const float &C1_STEPSIZE,
        const float &C2_STEPSIZE,
        const score::SaxsType::ScoreFunctionEnum &SCORE_FUNCTION,
        const bool &USE_ERRORS,
        const SaxsTransformation &TRANSFORMATION_TYPES,
        const bool &APPROXIMATE_SIDE_CHAINS,
        const bool &HARDWARE_TYPE
      );

      //! @brief Clone function
      //! @return pointer to new SaxsOptimization
      SaxsOptimization *Clone() const;

    /////////////////
    // data access //
    /////////////////

      const float &GetC1Min() const
      {
        return m_c1Min;
      }

      const float &GetC1Max() const
      {
        return m_c1Max;
      }

      const float &GetC2Min() const
      {
        return m_c2Min;
      }

      const float &GetC2Max() const
      {
        return m_c2Max;
      }

      const float &GetC1StepSize() const
      {
        return m_c1StepSize;
      }

      const float &GetC2StepSize() const
      {
        return m_c2StepSize;
      }

      void SetC1Min( const float &C1_MIN)
      {
        m_c1Min = C1_MIN;
      }

      void SetC1Max( const float &C1_MAX)
      {
        m_c1Max = C1_MAX;
      }

      void SetC2Min( const float &C2_MIN)
      {
        m_c2Min = C2_MIN;
      }

      void SetC2Max( const float &C2_MAX)
      {
        m_c2Max = C2_MAX;
      }

      void SetC1StepSize( const float &C1_STEPSIZE)
      {
        m_c1StepSize = C1_STEPSIZE;
      }

      void SetC2StepSize( const float &C2_STEPSIZE)
      {
        m_c2StepSize = C2_STEPSIZE;
      }

      void SetScoreFunction( const score::SaxsType::ScoreFunctionEnum &SCORE_FUNCTION)
      {
        m_scoreType = SCORE_FUNCTION;
      }

      void SetErrorFlag( const bool &USE_ERRORS)
      {
        m_useErrors = USE_ERRORS;
      }

      void SetTransformationTypes( const SaxsTransformation &TRANSFORMATION_TYPES)
      {
        m_transformation = TRANSFORMATION_TYPES;
      }

      void SetApproximateSideChainsFlag( const bool &APPROXIMATE_SIDE_CHAINS)
      {
        m_ApproximateSideChains = APPROXIMATE_SIDE_CHAINS;
      }

      void SetHardwareType( const bool &HARDWARE_TYPE)
      {
        m_Cpu = HARDWARE_TYPE;
      }

      //! @brief set the experimental saxs profile
      //! @param EXPERIMENTAL_DATA the new experimental data
      void SetExperimentalData( const util::ShPtr< SaxsScatteringData> &EXPERIMENTAL_DATA)
      {
        m_ExpData = EXPERIMENTAL_DATA;
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operators  //
    ////////////////

      //! @brief overloaded () operator to optimize c1 and c2 parameters
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return returns SaxsOptiResult: c1, c2, and min chi score
      SaxsOptiResult operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief returns lowest chi score from adusting c1 and c2 parameters over a grid of 3200 points
      //! @brief for each column of the grid each row is tested until the function at a given point increases
      //! @brief this method is faster than a complete grid search, but still slow to find a final answer
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return returns SaxsOptiResult: c1, c2, and min chi score
      SaxsOptiResult SmartGridSearch( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief returns lowest chi score from adusting c1 and c2 parameters over a grid of 3200 points
      //! @brief searches the grid by a line search using quadratic interpolation along each row of the grid
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return SaxsOptiResult: c1, c2, and min chi score
      SaxsOptiResult GridWalk( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief returns lowest chi score from a search pattern starting at quadratic interpolation line search
      //! @param RADIUS, search radius
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return SaxsOptiResult: c1, c2, and min chi score
      SaxsOptiResult HookJeeves
      (
        linal::VectorND< float, 2> INITIAL_POINT,
        float RADIUS,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief one Dimensional unimodal derivative free optimization method using GoldenSection Method
      //! @param MIN_C1 - Minimum C1 value
      //! @param MAX_C1 - Maximum C1 value
      //! @param C2_VALUE - Constant row of the matrix
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return returns storage vector of 3 values: c1, c2, and min chi score
      SaxsOptiResult GoldenSection
      (
        float MIN_C1,
        float MAX_C1,
        float C2_VALUE,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief one Dimensional unimodal derivative free optimization method using Quadratic Interpolation
      //! @param MIN_C1 - Minimum C1 value
      //! @param MAX_C1 - Maximum C1 value
      //! @param C2_VALUE - Constant row of the matrix
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return returns storage vector of 3 values: c1, c2, and min chi score
      SaxsOptiResult QuadraticInterpolation
      (
        float MIN_C1,
        float MAX_C1,
        float C2_VALUE,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    private:

      // helper functions

      //! @brief Convert C1 coordinates to Matrix position
      //! @param C1 - value of contrast density parameter
      //! @return matrix location of C1 value
      int C1ParamtoMatrix( float C1) const;

      //! @brief Convert C2 coordinates to Matrix position
      //! @param C2 - value of contrast density parameter
      //! @return matrix location of C2 value
      int C2ParamtoMatrix( float C2) const;

      //! @brief takes input point from line search method and rounds it to a valid start point on the grid
      //! @param INITIAL_POINT point from line search that is not on a grid point
      //! @return point on Grid
      linal::VectorND< float, 2> DiscretePoint( const linal::VectorND< float, 2> INITIAL_POINT)const;

      //! @brief compute chi for given parameters and protein model
      //! @param C1 - electron density adjustment parameter for form factors
      //! @param C2 - hydration shell adjustmet parameter for form factors
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return chi score between experimental and computed saxs profiles
      float EvaluateFunction( float C1, float C2, const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief returns lowest chi score all points around the center point based on Radius
      //! @param POINT, center point consisting of c1 and c2
      //! @param RADIUS, search radius
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return linal::VectorND< float, 2>: c1, c2 coordinates of min function value
      linal::VectorND< float, 2> Explore
      (
        const linal::VectorND< float, 2> &POINT,
        float RADIUS,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief returns the function value of a given c1 and c2 combination.  Results are stored in m_searchGrid to
      //! @brief prevent duplication of work
      //! @param POINT, center point consisting of c1 and c2
      //! @param PROTEIN_MODEL - protein model that contains the "experimental" data and coordinates for current model
      //! @return function value of given c1 and c2 parameter
      float EvaluatePoint( const linal::VectorND< float, 2> &POINT, const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief find index of minimum value in storage vector
      //! @param ARRAY storage vector to search
      //! @return index of minimum value in array
      size_t GetMinimumIndex( const storage::Vector< float> &ARRAY) const;

      //! @brief return minimum value in storage vector
      //! @param ARRAY storage vector to search
      //! @return minimum value in array
      float GetMinimumValue( const storage::Vector< float> &ARRAY) const;

      //! @brief Ensure that the proposed point to evaluate is inside the searchGrid
      //! @param POINT, center point consisting of c1 and c2
      //! @return true if inbounds
      bool Inbounds( linal::VectorND< float, 2> POINT) const;

      //! @brief Ensure that the proposed point to evaluate (D) is inside the specified interval
      //! @param A - Left Boundary Point
      //! @param B - Right Boundary Point
      //! @param TEST_POINT, point to evaluate
      //! @return true if inbounds
      bool InsideQuadraticInterval( float A, float B, float TEST_POINT) const;

      //! @brief computes a quadratic model based on 3 values and thier function
      //! @param A - First point
      //! @param FA - function at first point
      //! @param B - Second point
      //! @param FB - function at second point
      //! @param C - Third point
      //! @param FC - function at third point
      //! @return function point
      float Quadratic( float A, float FA, float B, float FB, float C, float FC) const;

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

    }; // class SaxsOptimization
  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAXS_OPTIMIZATION_H_
