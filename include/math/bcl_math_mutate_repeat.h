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

#ifndef BCL_MATH_MUTATE_REPEAT_H_
#define BCL_MATH_MUTATE_REPEAT_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_result.h"
#include "bcl_math_range.h"
#include "random/bcl_random_uniform_distribution.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateRepeat
    //! @brief template class MutateRepeat allows repetition of the same mutate on a given argument
    //! @details This class applies the given Mutate on the given argument N times consecutively, where N is specified
    //! by the user and returns the result. It makes sure that if a mutate leads to a undefined argument ( such as
    //! a skipped step in Monte Carlo minimization), then it continues iterations with the last defined argument.
    //!
    //! @see @link example_math_mutate_repeat.cpp @endlink
    //! @author woetzen
    //! @date Aug 23, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class MutateRepeat :
      public MutateInterface< t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! minimum number of repetitions
      size_t m_MinNumberRepeats;

      //! maximum number of repetitions
      size_t m_MaxNumberRepeats;

      //! ShPtr to Mutate to be repeated
      util::ShPtr< MutateInterface< t_ArgumentType> > m_Mutate;

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
      MutateRepeat< t_ArgumentType>() :
        m_MinNumberRepeats(),
        m_MaxNumberRepeats(),
        m_Mutate()
      {
      }

      //! @brief construct from a mutate and number of repeats
      //! @param SP_MUTATE ShPtr to the Mutate that is going to be repeated
      //! @param MIN_NUMBER_REPEATS minimum number of repetitions wanted
      MutateRepeat< t_ArgumentType>
      (
        const util::ShPtr< MutateInterface< t_ArgumentType> > &SP_MUTATE,
        const size_t MIN_NUMBER_REPEATS = 1,
        const size_t MAX_NUMBER_REPEATS = 0
      ) :
        m_MinNumberRepeats( MIN_NUMBER_REPEATS),
        m_MaxNumberRepeats( !MAX_NUMBER_REPEATS ? MIN_NUMBER_REPEATS : MAX_NUMBER_REPEATS),
        m_Mutate( SP_MUTATE)
      {
        BCL_Assert( m_MinNumberRepeats <= m_MaxNumberRepeats, "Min <= Max in MutateRepeat");
      }

      //! @brief Clone function
      //! @return pointer to new MutateRepeat< t_ArgumentType>
      MutateRepeat< t_ArgumentType> *Clone() const
      {
        return new MutateRepeat< t_ArgumentType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param ARGUMENT Argument that will undergo mutates
      //! @return MutateResult including the resultant argument
      MutateResult< t_ArgumentType> operator()( const t_ArgumentType &ARGUMENT) const
      {
        // make copy of argument that will be mutated
        util::ShPtr< t_ArgumentType> mutated_argument( ARGUMENT.Clone());

        const size_t desired_n_repeats
        (
          m_MinNumberRepeats != m_MaxNumberRepeats
          ? random::GetGlobalRandom().SizeT( Range< size_t>( m_MinNumberRepeats, m_MaxNumberRepeats + 1))
          : m_MinNumberRepeats
        );

        // iterate the number of specified repeats
        for( size_t repeat_no( 0); repeat_no != desired_n_repeats; ++repeat_no)
        {
          // mutate result
          MutateResult< t_ArgumentType> mutate_result( m_Mutate->operator()( *mutated_argument));

          // if the mutate was successful update the mutated argument
          if( mutate_result.GetArgument().IsDefined())
          {
            // update the mutated argument to the one return from this mutate behind mutate_itr
            mutated_argument = mutate_result.GetArgument();
          }
        }

        // instantiate the result that will be returned
        MutateResult< t_ArgumentType> result( mutated_argument, *this);

        // end
        return result;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_MinNumberRepeats, ISTREAM);
        io::Serialize::Read( m_MaxNumberRepeats, ISTREAM);
        io::Serialize::Read( m_Mutate, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_MinNumberRepeats, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_MaxNumberRepeats, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Mutate, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class MutateRepeat

    // instantiate s_Instance
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> MutateRepeat< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateRepeat< t_ArgumentType>())
    );

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_MUTATE_REPEAT_H_
