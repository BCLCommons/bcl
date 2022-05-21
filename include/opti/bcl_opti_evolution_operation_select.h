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

#ifndef BCL_OPTI_EVOLUTION_OPERATION_SELECT_H_
#define BCL_OPTI_EVOLUTION_OPERATION_SELECT_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_evolution.h"
#include "bcl_opti_evolution_population.h"
#include "bcl_opti_evolution_population_member.h"
#include "math/bcl_math_object_stochastic_selector.h"
#include "util/bcl_util_string_numeric_conversion.h"

namespace bcl
{
  namespace opti
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvolutionOperationSelect
    //! @brief probabilistically selects from a set of EvolutionOperations
    //!
    //! @see @link example_opti_evolution_operation_select.cpp @endlink
    //! @author geanesar
    //! @date Oct 15, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_MemberType, typename t_FitnessType>
      class EvolutionOperationSelect :
      public util::SerializableInterface
    {

    //////////////
    // Typedefs //
    //////////////
    public:

      //! Typedef for evolution operation functions
      typedef typename ApproximatorEvolution< t_MemberType, t_FitnessType>::EvolutionOperationSerializable EvolutionOperationSerializable;

    private:

      //! Decision node for picking operations
      math::ObjectStochasticSelector< EvolutionOperationSerializable> m_DecisionNode;

    public:

      //! A single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief Default constructor
      EvolutionOperationSelect() :
        m_DecisionNode()
      {
      }

      //! @brief select an operation
      //! @return a reference to the operation
      const EvolutionOperationSerializable &Select() const
      {
        return m_DecisionNode.SelectRandomCase();
      }

      //! @brief Clone function
      //! @return new Pointer to a copy of the actual object behind the pointer
      EvolutionOperationSelect *Clone() const
      {
        return new EvolutionOperationSelect( *this);
      }

      //! @brief returns class name of this object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of this object when used in a dynamic context
      //! @return the name of this object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "OperationSelect");
        return s_alias;
      }

      //! @brief add an operation to the decision selector
      void AddImplementation( const EvolutionOperationSerializable &OPERATION, const double &PROBABILITY)
      {
        m_DecisionNode.AddImplementation( OPERATION, PROBABILITY);
      }

      //! @brief returns the decision node used by this class
      //! @return the decision node
      math::ObjectStochasticSelector< EvolutionOperationSerializable> &GetDecisionNode()
      {
        return m_DecisionNode;
      }

      //! @brief get the size of the decision node
      //! @return the size of the decision node
      size_t GetSize() const
      {
        return m_DecisionNode.GetSize();
      }

      //! @brief called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
      {
        return m_DecisionNode.ReadInitializerSuccessHook( SERIALIZER, ERR_STREAM);
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {

        io::Serializer member_data;

        member_data.AddInitializer
        (
          "",
          "Mutation operations and probabilities; specify them as: (objects=(...),probabilities=(...))",
          io::Serialization::GetAgent( &m_DecisionNode)
        );

        member_data.SetClassDescription( "Select and use a probabilistically-weighted random operation from a set of operations");

        return member_data;
      }

    };
  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_EVOLUTION_OPERATION_SELECT_H_
