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

#ifndef BCL_MATH_MUTATE_DECISION_NODE_H_
#define BCL_MATH_MUTATE_DECISION_NODE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_interface.h"
#include "bcl_math_mutate_result.h"
#include "bcl_math_object_probability_distribution.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateDecisionNode
    //! @brief MutateDecisionNode is an adapter for ObjectProbabilityDistribution to be used in MutateInterface framework
    //! @details This MutateInterface derived class, stores internally an object probability distribution that contains
    //! set of mutates and associated probabilities.This class forms the root and the other nodes in the mutate tree.
    //!
    //! @see @link example_math_mutate_decision_node.cpp @endlink
    //! @author karakam, woetzen
    //! @date Nov 20, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class MutateDecisionNode :
      public MutateInterface< t_DataType>
    {

    //////////
    // data //
    //////////

    private:

      //! list of function interfaces for mutating with their corresponding probabilities
      ObjectProbabilityDistribution< MutateInterface< t_DataType> > m_Distribution;

      //! string alias
      std::string m_Alias;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateDecisionNode() :
        m_Distribution()
      {
      }

      //! @brief construct from ObjectProbabilityDistribution and a scheme
      //! @param DISTRIBUTION ObjectProbabilityDistribution
      //! @param SCHEME scheme to be used
      explicit MutateDecisionNode
      (
        const ObjectProbabilityDistribution< MutateInterface< t_DataType> > &DISTRIBUTION,
        const std::string &ALIAS = std::string( "MutateDecisionNode")
      ) :
        m_Distribution( DISTRIBUTION),
        m_Alias( ALIAS)
      {
      }

      //! @brief virtual copy constructor
      MutateDecisionNode< t_DataType> *Clone() const
      {
        return new MutateDecisionNode< t_DataType>( *this);
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

      //! @brief alias
      //! @return an alias for the class
      const std::string &GetAlias() const
      {
        return m_Alias;
      }

      //! @brief object probability distribution
      //! @return object probability distribution
      const ObjectProbabilityDistribution< MutateInterface< t_DataType> > &GetObjectProbabilityDistribution() const
      {
        return m_Distribution;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer parameters;
        parameters.SetClassDescription( "Stores an object probability distribution.");
        parameters.AddInitializer
        (
          "probability distribution",
          "probability distribution of the objects",
          io::Serialization::GetAgent( &m_Distribution)
        );

        return parameters;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add muate with probability
      //! @param MUTATE the mutate to add
      //! @param PROBABILITY the probability for that mutate
      void AddMutate( const MutateInterface< t_DataType> &MUTATE, const double PROBABILITY)
      {
        m_Distribution.PushBack( PROBABILITY, MUTATE);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief Calls the operator of a randomly determine Mutate and returns the result for given ARGUMENT
      //! @param ARGUMENT Argument of interest
      MutateResult< t_DataType> operator()( const t_DataType &ARGUMENT) const
      {
        // determine random case
        util::SiPtr< const MutateInterface< t_DataType> > selected_mutate( m_Distribution.DetermineRandomCase());

        // issue message
        BCL_MessageVrb( "Selected: " + selected_mutate->GetScheme());

        // get mutate result
        MutateResult< t_DataType> result( selected_mutate->operator()( ARGUMENT));

        // insert this node in the mutate node list
        result.AddNode( *this);

        // return result
        return result;
      }

    //////////////////////
    // input and output //
    //////////////////////

    }; // template class MutateDecisionNode

    //! single instance of this class
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> MutateDecisionNode< t_DataType>::s_Instance
    (
      util::Enumerated< MutateInterface< t_DataType> >::AddInstance( new MutateDecisionNode< t_DataType>())
    );

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_MUTATE_DECISION_NODE_H_
