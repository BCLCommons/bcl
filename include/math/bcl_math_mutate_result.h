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

#ifndef BCL_MATH_MUTATE_RESULT_H_
#define BCL_MATH_MUTATE_RESULT_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateResult
    //! @brief provides functionality to return mutated argument as well as Mutates that were called
    //! @details The MutateResult class, stores a pointer to the mutated object, so that if the mutate function fails,
    //! there is no need to make a copy of the unchanged argument. It also provides list of MutateInterface that have
    //! been called in reverse order so that the call can be traced back if necessary. This SiPtrList of MutateInterface
    //! is expected to be filled as the result travels back from the leaf node in the tree structure back to the root
    //! thus back to the Approximater of interest. This bookkeeping allows tracking the success rates of individual
    //! Mutate nodes as well as paths to such nodes which can be used for various functionalities such as success based
    //! weight adjusting for Mutate nodes
    //!
    //! @tparam t_ArgumentType the type of the argument
    //!
    //! @see @link example_math_mutate_result.cpp @endlink
    //! @author karakam, woetzen
    //! @date 08.15.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class MutateResult :
      public util::SerializableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! ShPtr to the mutated argument
      util::ShPtr< t_ArgumentType> m_Argument;

      //! List of Mutates that describes the path traveled from Mutate node from which this MutateResult originated
      util::SiPtrList< const MutateInterface< t_ArgumentType> > m_Nodes;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateResult() :
        m_Argument(),
        m_Nodes()
      {
      }

      //! @brief constructor from an argument and a mutate node of origin
      //! @param ARGUMENT ShPtr to mutated argument to be returned as result
      //! @param NODE Mutate node which produced the result
      MutateResult
      (
        const util::ShPtr< t_ArgumentType> &ARGUMENT,
        const MutateInterface< t_ArgumentType> &NODE
      ) :
        m_Argument( ARGUMENT),
        m_Nodes( 1, &NODE)
      {
      }

      //! @brief constructor from an argument and a list of mutates
      //! @param ARGUMENT ShPtr to mutated argument to be returned as result
      //! @param MUTATE_LIST list of mutates that produce the result
      MutateResult
      (
        const util::ShPtr< t_ArgumentType> &ARGUMENT,
        const util::SiPtrList< const MutateInterface< t_ArgumentType> > &MUTATE_LIST
      ) :
        m_Argument( ARGUMENT),
        m_Nodes( MUTATE_LIST)
      {
      }

      //! @brief virtual copy constructor
      MutateResult< t_ArgumentType> *Clone() const
      {
        return new MutateResult< t_ArgumentType>( *this);
      };

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "Iterations");
        return s_alias;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
         io::Serializer serializer;
         serializer.SetClassDescription( "Triggers after the given number of iterations");
         return serializer;
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return mutated argument
      //! @return mutated argument
      const util::ShPtr< t_ArgumentType> &GetArgument() const
      {
        return m_Argument;
      }

      //! @brief return mutate nodes
      //! @return mutate nodes
      const util::SiPtrList< const MutateInterface< t_ArgumentType> > &GetNodes() const
      {
        return m_Nodes;
      }

      //! @brief insert a MutateInterface into the nodes at the front
      //! @param MUTATE_INTERFACE mutate that was contributing to the argument result
      void AddNode( const MutateInterface< t_ArgumentType> &MUTATE_INTERFACE)
      {
        m_Nodes.PushFront( util::ToSiPtr( MUTATE_INTERFACE));
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator += for adding a mutate result to this
      //! @param MUTATE_RESULT the MutateResult that will be added to this
      //! @return copy of this
      MutateResult< t_ArgumentType> operator +=( const MutateResult< t_ArgumentType> &MUTATE_RESULT)
      {
        // return this
        return *this;
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
        ISTREAM >> m_Argument
                >> m_Nodes;

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
        io::Serialize::Write( m_Argument, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Nodes, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class MutateResult

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_MUTATE_RESULT_H_
