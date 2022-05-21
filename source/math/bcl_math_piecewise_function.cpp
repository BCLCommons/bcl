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
#include "math/bcl_math_piecewise_function.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PiecewiseFunction::s_Instance
    (
      GetObjectInstances().AddInstance( new PiecewiseFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PiecewiseFunction::PiecewiseFunction() :
      m_RangesFunctions()
    {
    }

    //! @brief construct from the ranges and their functions the piecewise function
    //! @param RANGES_FUNCTIONS List which gives the ranges and their functions as Pairs
    PiecewiseFunction::PiecewiseFunction
    (
      const storage::List
      <
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
      > &RANGES_FUNCTIONS
    ) :
      m_RangesFunctions( RANGES_FUNCTIONS)
    {
      // check to make sure the ranges and functions define a valid piecewise function (i.e. no overlapping ranges
      // or no non existing functions)
      BCL_Assert
      (
        RangesFunctionsValid( m_RangesFunctions),
        "The ShPtr provided does not point to anything or the ranges overlap."
      );
    }

    //! @brief Clone function
    //! @return pointer to new class_name
    PiecewiseFunction *PiecewiseFunction::Clone() const
    {
      return new PiecewiseFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PiecewiseFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetRangesFunctions returns m_RangesFunctions
    //! @return returns List which is m_RangesFunctions
    const storage::List
    <
      storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
    > &PiecewiseFunction::GetRangesFunctions() const
    {
      return m_RangesFunctions;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking an x-value and returning the y-value based on "m_RangesFunctions"
    //! @param X_VALUE double which is the x value from which the y-value will be calculated
    //! @return returns a double which is the y-value based on X_VALUE and "m_RangesFunctions"
     double PiecewiseFunction::operator()( const double &X_VALUE) const
     {
       // iterate through the list of ranges and functions to see which range the X_VALUE lies within
       for
       (
         storage::List
         <
           storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
         >::const_iterator iter( m_RangesFunctions.Begin()), iter_end( m_RangesFunctions.End());
         iter != iter_end; ++iter
       )
       {
         const Range< double> &range_a( iter->First());
         const util::ShPtr< FunctionInterfaceSerializable< double, double> > &function( iter->Second());

         // true if X_VALUE is within the range currently denoted by "iter"
         if( range_a.IsWithin( X_VALUE))
         {
           // execute the function to get the corresponding function value
           return function->operator()( X_VALUE);
         }
       }

       // if this point is reached then "X_VALUE" does not fall in any of the ranges so return undefined double
       return util::GetUndefined< double>();
     }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PiecewiseFunction::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RangesFunctions, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PiecewiseFunction::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RangesFunctions, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determine whether the ranges overlap or if the ShPtr of the functions point to anything
    //! @param RANGES_FUNCTIONS_VALID gives a List with a pair of ranges and function interfaces to be analyzed
    //! @return bool which will tell whether or not the list is valid
    bool PiecewiseFunction::RangesFunctionsValid
    (
      const storage::List
      <
        storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
      > &RANGES_FUNCTIONS_VALID
    )
    {
      // iterate through the list "RANGES_FUNCTIONS_VALID"
      for
      (
        storage::List
        <
          storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
        >::const_iterator iter_a( RANGES_FUNCTIONS_VALID.Begin()), iter_end( RANGES_FUNCTIONS_VALID.End());
        iter_a != iter_end; ++iter_a
      )
      {
        const Range< double> &range_a( iter_a->First());
        const util::ShPtr< FunctionInterfaceSerializable< double, double> > &function( iter_a->Second());

        // for each function determine if it exists
        if( !function.IsDefined())
        {
          BCL_MessageCrt( "The function does not exist.")
          return false;
        }

        // iterate through the list "RANGES_FUNCTIONS_VALID" again and
        // for each range determine if it overlaps with any of the other ranges
        for
        (
          storage::List
          <
            storage::Pair< Range< double>, util::ShPtr< FunctionInterfaceSerializable< double, double> > >
          >::const_iterator iter_b( iter_a); iter_b != iter_end; ++iter_b
        )
        {
          const Range< double> &range_b( iter_b->First());

          // ignore the ranges that are equal to themselves
          if( iter_b == iter_a)
          {
            continue;
          }

          // true if "range_a" overlaps with "range_b" or vice versa
          if( range_a.DoesOverlap( range_b))
          {
            BCL_MessageCrt
            (
              "The following ranges overlap: "
                + range_a.GetString( util::Format()) + ", " + range_b.GetString( util::Format()) + "."
            );

            // return false indicating that the list of ranges and functions is not valid
            return false;
          }
        }
      }

      // return true indicating that none of the ranges overlap and all of the functions are defined
      return true;
    }

  } // namespace math
} // namespace bcl
