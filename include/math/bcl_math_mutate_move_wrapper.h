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

#ifndef BCL_MATH_MUTATE_MOVE_WRAPPER_H_
#define BCL_MATH_MUTATE_MOVE_WRAPPER_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_interface.h"
#include "coord/bcl_coord_move_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateMoveWrapper
    //! @brief wraps a move interface in a MutateMoveWrapper derived class
    //! @details This class allows using a move directly in place of MutateMoveWrapper derived class
    //!
    //! @see @link example_math_mutate_move_wrapper.cpp @endlink
    //! @author karakam
    //! @date Mar 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_ArgumentType>
    class MutateMoveWrapper :
      public MutateInterface< t_ArgumentType>
    {
    private:

    //////////
    // data //
    //////////

      //! ShPtr to the MoveInterface derived class that is wrapper
      util::Implementation< coord::MoveInterface> m_Move;

      //! whether to do a hardcopy( just a clone otherwise)
      bool m_DoHardCopy;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateMoveWrapper() :
        m_Move(),
        m_DoHardCopy( true)
      {
      }

      //! @brief constructor from a MoveInterface
      //! @param MOVE MoveInterface derived class to be wrapped
      //! @param DO_HARD_COPY  whether to do a hardcopy( just a clone otherwise)
      MutateMoveWrapper
      (
        const coord::MoveInterface &MOVE,
        const bool DO_HARD_COPY = true
      ) :
        m_Move( MOVE),
        m_DoHardCopy( DO_HARD_COPY)
      {
      }

      //! @brief constructor from a ShPtr to a MoveInterface
      //! @param SP_MOVE ShPtr to MoveInterface derived class to be wrapped
      //! @param DO_HARD_COPY  whether to do a hardcopy( just a clone otherwise)
      MutateMoveWrapper
      (
        const util::ShPtr< coord::MoveInterface> &SP_MOVE,
        const bool DO_HARD_COPY = true
      ) :
        m_Move( *SP_MOVE),
        m_DoHardCopy( DO_HARD_COPY)
      {
      }

      //! @brief virtual copy constructor
      //! @return pointer to a new MutateMoveWrapper copied from this one
      MutateMoveWrapper< t_ArgumentType> *Clone() const
      {
        return new MutateMoveWrapper< t_ArgumentType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "MutateMoveWrapper");
        return s_name;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Applies a random rotation relative to an external coordinate system.");
        serializer.AddInitializer
        (
          "move interface",
          "move interface that is wrapped in this class",
          io::Serialization::GetAgent( &m_Move)
        );
        serializer.AddInitializer
        (
          "hardcopy",
          "make a hard copy as opposed to cloning",
          io::Serialization::GetAgent( &m_DoHardCopy)
        );

        return serializer;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param ARGUMENT Argument of interest
      //! @return MutateResult that results from mutating to the argument
      MutateResult< t_ArgumentType> operator()( const t_ArgumentType &ARGUMENT) const
      {
        // make a clone or hardcopy of the argument depending on member variable m_DoHardCopy
        util::ShPtr< t_ArgumentType> sp_argument
        (
          m_DoHardCopy ? ARGUMENT.HardCopy() : ARGUMENT.Clone()
        );

        // apply the move
        m_Move->Move( *sp_argument);

        // apply the move to the argument and return the result in a ShPtr
        return MutateResult< t_ArgumentType>( sp_argument, *this);
      }

    //////////////////////
    // input and output //
    //////////////////////

    }; // template class MutateMoveWrapper

    // instantiate s_Instance
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> MutateMoveWrapper< t_ArgumentType>::s_Instance
    (
      util::Enumerated< MutateInterface< t_ArgumentType> >::AddInstance( new MutateMoveWrapper< t_ArgumentType>())
    );

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_MUTATE_MOVE_WRAPPER_H_
