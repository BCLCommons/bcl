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

#ifndef BCL_MATH_MUTATE_COMBINE_H_
#define BCL_MATH_MUTATE_COMBINE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_result.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateCombine
    //! @brief template class MutateCombine applied multiple mutates to the given argument consecutively
    //! @details This mutate class, initialized with ShPtrList of mutates, applies each of the given mutates
    //! consecutively to the given argument.
    //!
    //! @see @link example_math_mutate_combine.cpp @endlink
    //! @author woetzen
    //! @date Aug 23, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class MutateCombine :
      public MutateInterface< t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! list of mutate functions
      util::ShPtrList< FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > > m_Mutates;

      //! if true the mutates will continue with the last defined result
      bool m_IgnoreUndefinedResult;

      //! scheme of this mutate
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme()
      {
        static const std::string s_scheme( "mutate_combine");
        return s_scheme;
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateCombine( const std::string &SCHEME = GetDefaultScheme()) :
        m_Mutates(),
        m_IgnoreUndefinedResult(),
        m_Scheme( SCHEME)
      {
      }

      //! @brief construct from list of MutateInterfaces
      //! @param MUTATES list of MutateInterfaces that will be use to make member variable
      MutateCombine
      (
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_A,
        bool IGNORE_UNDEFINED_RESULT,
        const std::string &SCHEME
      ) :
        m_Mutates(),
        m_IgnoreUndefinedResult( IGNORE_UNDEFINED_RESULT),
        m_Scheme( SCHEME)
      {
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_A));
      }

      //! @brief construct from list of MutateInterfaces
      //! @param MUTATES list of MutateInterfaces that will be use to make member variable
      MutateCombine
      (
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_A,
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_B,
        bool IGNORE_UNDEFINED_RESULT,
        const std::string &SCHEME
      ) :
        m_Mutates(),
        m_IgnoreUndefinedResult( IGNORE_UNDEFINED_RESULT),
        m_Scheme( SCHEME)
      {
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_A));
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_B));
      }

      //! @brief construct from list of MutateInterfaces
      //! @param MUTATES list of MutateInterfaces that will be use to make member variable
      MutateCombine
      (
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_A,
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_B,
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_C,
        bool IGNORE_UNDEFINED_RESULT,
        const std::string &SCHEME
      ) :
        m_Mutates(),
        m_IgnoreUndefinedResult( IGNORE_UNDEFINED_RESULT),
        m_Scheme( SCHEME)
      {
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_A));
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_B));
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_C));
      }

      //! @brief construct from list of MutateInterfaces
      //! @param MUTATES list of MutateInterfaces that will be use to make member variable
      MutateCombine
      (
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_A,
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_B,
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_C,
        const FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > &MUTATE_D,
        bool IGNORE_UNDEFINED_RESULT,
        const std::string &SCHEME
      ) :
        m_Mutates(),
        m_IgnoreUndefinedResult( IGNORE_UNDEFINED_RESULT),
        m_Scheme( SCHEME)
      {
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_A));
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_B));
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_C));
        m_Mutates.PushBack( util::CloneToShPtr( MUTATE_D));
      }

      //! @brief construct from list of MutateInterfaces
      //! @param MUTATES list of MutateInterfaces that will be use to make member variable
      MutateCombine
      (
        const util::ShPtrList< MutateInterface< t_ArgumentType> > &MUTATES,
        bool IGNORE_UNDEFINED_RESULT,
        const std::string &SCHEME
      ) :
        m_Mutates(),
        m_IgnoreUndefinedResult( IGNORE_UNDEFINED_RESULT),
        m_Scheme( SCHEME)
      {
        // iterate through the mutates to fill m_Mutates
        for
        (
          typename util::ShPtrList< MutateInterface< t_ArgumentType> >::const_iterator
            itr( MUTATES.Begin()), itr_end( MUTATES.End());
          itr != itr_end; ++itr
        )
        {
          // add current mutate to m_Mutates
          m_Mutates.PushBack( *itr);
        }
      }

      //! @brief construct from list of function interfaces
      //! @param MUTATES list of function interfaces that will be use to make member variable
      MutateCombine
      (
        const util::ShPtrList< FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> > > &MUTATES,
        bool IGNORE_UNDEFINED_RESULT,
        const std::string &SCHEME
      ) :
        m_Mutates( MUTATES),
        m_IgnoreUndefinedResult( IGNORE_UNDEFINED_RESULT),
        m_Scheme( SCHEME)
      {
      }

      //! @brief Clone function
      //! @return pointer to new MutateCombine< t_ArgumentType>
      MutateCombine< t_ArgumentType> *Clone() const
      {
        return new MutateCombine< t_ArgumentType>( *this);
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

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
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

        bool had_defined( false);
        util::ShPtr< t_ArgumentType> last_defined_mutated_argument( mutated_argument);

        // apply all mutations to the mutated argument
        for
        (
          typename util::ShPtrList
          <
            FunctionInterfaceSerializable< t_ArgumentType, MutateResult< t_ArgumentType> >
          >::const_iterator mutate_itr( m_Mutates.Begin()), mutate_itr_end( m_Mutates.End());
          mutate_itr != mutate_itr_end && ( mutated_argument.IsDefined() || m_IgnoreUndefinedResult);
          ++mutate_itr
        )
        {
          // update the mutated argument to the one return from this mutate behind mutate_itr
          mutated_argument = ( ( *mutate_itr)->operator()( *last_defined_mutated_argument)).GetArgument();

          if( mutated_argument.IsDefined())
          {
            had_defined = true;
            last_defined_mutated_argument = mutated_argument;
          }
        }

        // end
        return MutateResult< t_ArgumentType>
               (
                 had_defined || m_Mutates.IsEmpty() ? last_defined_mutated_argument : util::ShPtr< t_ArgumentType>(),
                 *this
               );
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
        // read member
        io::Serialize::Read( m_Mutates, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_Mutates, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class MutateCombine

    // instantiate s_Instance
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> MutateCombine< t_ArgumentType>::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateCombine< t_ArgumentType>())
    );

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_MUTATE_COMBINE_H_
