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

#ifndef BCL_OPTI_ENSEMBLE_NODE_H_
#define BCL_OPTI_ENSEMBLE_NODE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_optimization_interface.h"
#include "assemble/bcl_assemble_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnsembleNode
    //! @brief Optimization of data ensembles through optimization of its data members
    //!
    //! @see @link example_opti_ensemble_node.cpp @endlink
    //! @author fischea
    //! @date Oct 30, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class EnsembleNode :
      public OptimizationInterface< assemble::Ensemble< t_DataType> >
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! optimization implementation to apply on the members of the ensemble
      util::Implementation< OptimizationInterface< t_DataType> > m_Optimizer;

      //! number of independent optimization trajectories per ensemble member
      size_t m_NumberTrajectories;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief default constructor
      EnsembleNode() :
        m_Optimizer(),
        m_NumberTrajectories()
      {
      }

      //! @brief construct from optimization implementation
      //! @param OPTIMIZER algorithm that will be applied to each member of the ensemble
      //! @param NUMBER_TRAJECTORIES number of independent optimization trajectories per ensemble member
      EnsembleNode( const OptimizationInterface< t_DataType> &OPTIMIZER, size_t NUMBER_TRAJECTORIES = 1) :
        m_Optimizer( OPTIMIZER),
        m_NumberTrajectories( NUMBER_TRAJECTORIES)
      {
      }

      //! @brief clone function
      //! @return pointer to a new EnsembleNode
      EnsembleNode *Clone() const
      {
        return new EnsembleNode( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "EnsembleNode");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Ensemble Node.");
        serializer.AddInitializer
        (
          "optimizer",
          "algorithm applied to each member of the ensemble",
          io::Serialization::GetAgent( &m_Optimizer)
        );
        serializer.AddInitializer
        (
          "number trajectories",
          "number of independent optimization trajectories per ensemble member",
          io::Serialization::GetAgent( &m_NumberTrajectories),
          "1"
        );

        return serializer;
      }

    ///////////////
    // operations//
    ///////////////

    protected:

      //! @brief performs pre-processing of the argument before the optimization is performed
      //! @detail increases the ensemble by copying models if multiple models should be sampled per ensemble member
      //! @param ARGUMENT data to be pre-processed
      //! @return pre-processed data
      void PreProcess( assemble::Ensemble< t_DataType> &ARGUMENT) const
      {
        if( m_NumberTrajectories > 1)
        {
          storage::Vector< util::ShPtr< t_DataType> > new_members;
          for( auto member_it( ARGUMENT.Begin()), member_it_end( ARGUMENT.End()); member_it != member_it_end; ++member_it)
          {
            for( size_t copy_number( 1); copy_number < m_NumberTrajectories; ++copy_number)
            {
              util::ShPtr< t_DataType> member_copy( member_it->GetElement().HardCopy());
              new_members.PushBack( member_copy);
            }
          }

          // add the new members to the ensemble
          for( auto member_it( new_members.Begin()), member_it_end( new_members.End()); member_it != member_it_end; ++member_it)
          {
            ARGUMENT.AddElement( **member_it);
          }
        }
      }

      //! @brief initiates the optimization
      //! @param ARGUMENT data to be optimized
      void Optimize( assemble::Ensemble< t_DataType> &ARGUMENT) const
      {
        for( auto member_it( ARGUMENT.Begin()), member_it_end( ARGUMENT.End()); member_it != member_it_end; ++member_it)
        {
          ( *m_Optimizer)( member_it->GetElement());
        }
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class EnsembleNode

    //! single instance of this class
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> EnsembleNode< t_ArgumentType>::s_Instance
    (
      util::Enumerated< OptimizationInterface< assemble::Ensemble< t_ArgumentType> > >::AddInstance
      (
        new EnsembleNode< t_ArgumentType>()
      )
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_ENSEMBLE_NODE_H_
