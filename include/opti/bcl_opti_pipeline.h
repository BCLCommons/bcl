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

#ifndef BCL_OPTI_PIPELINE_H_
#define BCL_OPTI_PIPELINE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_optimization_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Pipeline
    //! @brief Pipeline for data optimization
    //! @detail
    //!
    //! @see @link example_opti_pipeline.cpp @endlink
    //! @author fischea
    //! @date Oct 24, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_ArgumentType>
    class Pipeline :
      public OptimizationInterface< t_ArgumentType>
    {

    //////////////
    // typedefs //
    //////////////

      //! module of the pipeline
      typedef util::Implementation< OptimizationInterface< t_ArgumentType> > PipelineModule;

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! modules forming the pipeline
      storage::Vector< PipelineModule> m_Pipeline;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief default constructor
      Pipeline() :
        m_Pipeline()
      {
      }

      //! @brief construct from data members
      Pipeline( const storage::Vector< PipelineModule> &PIPELINE) :
        m_Pipeline( PIPELINE)
      {
      }

      //! @brief clone function
      //! @return pointer to a new Pipeline
      Pipeline *Clone() const
      {
        return new Pipeline( *this);
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
        static const std::string s_alias( "Pipeline");
        return s_alias;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Pipeline for data optimization.");
        serializer.AddInitializer
        (
          "pipeline",
          "elements of the pipeline",
          io::Serialization::GetAgent( &m_Pipeline)
        );

        return serializer;
      }

      //! @brief Append a module to the pipeline
      //! @param MODULE module to append to the pipeline
      void AppendModule( const OptimizationInterface< t_ArgumentType> &MODULE)
      {
        m_Pipeline.PushBack( MODULE);
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief applies the pipeline on the given argument
      //! @param ARGUMENT data on which the pipeline will be applied
      void Optimize( t_ArgumentType &ARGUMENT) const
      {
        // iterate over the pipeline and apply each optimizer to the argument
        for( auto it_module( m_Pipeline.Begin()), it_module_end( m_Pipeline.End()); it_module != it_module_end; ++it_module)
        {
          ( **it_module)( ARGUMENT);
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class Pipeline

    //! single instance of this class
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> Pipeline< t_ArgumentType>::s_Instance
    (
      util::Enumerated< OptimizationInterface< t_ArgumentType> >::AddInstance( new Pipeline< t_ArgumentType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_PIPELINE_H_
