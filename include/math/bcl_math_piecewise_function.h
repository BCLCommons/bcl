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

#ifndef BCL_MATH_PIECEWISE_FUNCTION_H_
#define BCL_MATH_PIECEWISE_FUNCTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "bcl_math_range.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PiecewiseFunction
    //! @brief calculates piecewise functions from a given "X" value
    //! @details http://en.wikipedia.org/wiki/Piecewise
    //!
    //! @see @link example_math_piecewise_function.cpp @endlink
    //! @author akinlr, alexanns
    //! @date Jun 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PiecewiseFunction :
      public FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      //! Stores the ranges and corresponding functions needed for the piecewise function
      storage::List
      <
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
      > m_RangesFunctions;

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
      PiecewiseFunction();

      //! @brief construct from the ranges and their functions the piecewise function
      //! @param RANGES_FUNCTIONS List which gives the ranges and their functions as Pairs
      PiecewiseFunction
      (
        const storage::List
        <
          storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        > &RANGES_FUNCTIONS
      );

      //! @brief Clone function
      //! @return pointer to new class_name
      PiecewiseFunction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetRangesFunctions return m_RangesFunctions
      //! @return returns List which is m_RangesFunctions
      const storage::List
      <
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
      > &GetRangesFunctions() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() taking an x-value and returning the y-value based on "m_RangesFunctions"
      //! @param X_VALUE double which is the x-value from which the y-value will be calculated
      //! @return returns a double which is the y-value based on X_VALUE and "m_RangesFunctions"
      double operator()( const double &X_VALUE) const;

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

      //! @brief determine whether the ranges overlap or if the ShPtr of the functions point to anything
      //! @param RANGES_FUNCTIONS_VALID gives a List with a pair of ranges and function interfaces to be analyzed
      //! @return bool which will tell whether or not the list is valid
      static bool RangesFunctionsValid
      (
        const storage::List
        <
          storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        > &RANGES_FUNCTIONS_VALID
      );

    }; // class PiecewiseFunction

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_PIECEWISE_FUNCTION_H_
