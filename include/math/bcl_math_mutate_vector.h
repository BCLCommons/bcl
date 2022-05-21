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

#ifndef BCL_MATH_MUTATE_VECTOR_H_
#define BCL_MATH_MUTATE_VECTOR_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_interface.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateVector
    //! @brief mutates a vector of parameters
    //! @details allows you to give an acceptable range for mutation of each parameter in a given vector, a random
    //! number within this range will be added/subtracted to each parameter. The ranges can also be specifed as
    //! percentages instead of absolute values
    //!
    //! @see @link example_math_mutate_vector.cpp @endlink
    //! @author woetzen, durhamea, karakam
    //! @date 10/30/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateVector :
      public MutateInterface< linal::Vector< double> >
    {

    private:

    //////////
    // data //
    //////////

      //! the maximal change for each parameter during this mutation
      linal::Vector< double> m_Range;

      //! the number of params mutated at each call to the operator function
      size_t m_NumberParamsMutated;

      //! boolean to indicate whether the mutate ranges are percentages or not
      bool m_UsePercentage;

      //! boolean wether to prevent that the numbers in the vector become negative
      bool m_KeepPositive;

      //! index in the vector to be kept constant. If undefined, no index is protected.
      //! having such an index should prevent trivially different vectors due to scaling
      size_t m_Constant;

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
      MutateVector();

      //! @brief constructor from a fixed range, total number params, number params mutated and a boolean for using percentages
      //! @param NUMBER_PARAMS Total number of parameters
      //! @param MUTATE_RANGE the range allowed for all parameters
      //! @param NUMBER_PARAMS_MUTATED the number of params mutated at each call to the operator function
      //! @param USE_PERCENTAGE boolean to indicate whether the mutate ranges are percentages
      //! @param KEEP_POSITIVE keep the values as positive numbers
      //! @param CONSTANT index to keep constant. By default, no such index is defined
      MutateVector
      (
        const size_t NUMBER_PARAMS,
        const double MUTATE_RANGE,
        const size_t NUMBER_PARAMS_MUTATED,
        const bool USE_PERCENTAGE,
        const bool KEEP_POSITIVE,
        const size_t &CONSTANT = util::GetUndefinedSize_t()
      );

      //! @brief constructor from mutate range, number params mutated and a boolean for using percentages
      //! @param MUTATE_RANGES a vector of acceptable ranges for each parameter
      //! @param NUMBER_PARAMS_MUTATED the number of params mutated at each call to the operator function
      //! @param USE_PERCENTAGE boolean to indicate whether the mutate ranges are percentages
      //! @param KEEP_POSITIVE keep the values as positive numbers
      //! @param CONSTANT index to keep constant. By default, no such index is defined
      MutateVector
      (
        const linal::Vector< double> &MUTATE_RANGES,
        const size_t NUMBER_PARAMS_MUTATED,
        const bool USE_PERCENTAGE,
        const bool KEEP_POSITIVE,
        const size_t &CONSTANT = util::GetUndefinedSize_t()
      );

      //! @brief Clone function
      //! @return pointer to new MutateVector
      MutateVector *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an VECTOR and returning a mutated Vector
      //! @param VECTOR Vector of parameters to be mutated
      //! @return MutateResult containing the mutated Vector
      MutateResult< linal::Vector< double> > operator()( const linal::Vector< double> &VECTOR) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // MutateVector

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_MUTATE_VECTOR_H_
