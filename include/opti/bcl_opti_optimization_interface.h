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

#ifndef BCL_OPTI_OPTIMIZATION_INTERFACE_H_
#define BCL_OPTI_OPTIMIZATION_INTERFACE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_processor_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OptimizationInterface
    //! @brief Interface for optimization methods. This interface was introduced to provide an abstraction layer
    //! for optimization methods.
    //!
    //! @remarks example unnecessary
    //! @author fischea
    //! @date Aug 28, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_ArgumentType>
    class OptimizationInterface :
      public util::SerializableInterface
    {

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief virtual copy constructor
      //! @return pointer to a new OptimizationInterface
      virtual OptimizationInterface *Clone() const = 0;

    private:

      //! processors applied to the optimization target before initiating the optimization
      storage::Vector< util::Implementation< ProcessorInterface< t_ArgumentType> > > m_PreProcessors;

      //! processors applied to the optimization target after conclusion of the optimization
      storage::Vector< util::Implementation< ProcessorInterface< t_ArgumentType> > > m_PostProcessors;

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief add the provided algorithm to the pre-processors
      //! @param PROCESSOR algorithm to add to the pre-processors
      void AddPreProcessor( const ProcessorInterface< t_ArgumentType> &PROCESSOR)
      {
        m_PreProcessors.PushBack( PROCESSOR);
      }

      //! @brief add the provided algorithm to the post-processors
      //! @param PROCESSOR algorithm to add to the post-processors
      void AddPostProcessor( const ProcessorInterface< t_ArgumentType> &PROCESSOR)
      {
        m_PostProcessors.PushBack( PROCESSOR);
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "OptimizationInterface");
        serializer.AddInitializer
        (
          "preprocessors",
          "algorithms applied to the argument before the optimization is started",
          io::Serialization::GetAgent( &m_PreProcessors),
          "()"
        );
        serializer.AddInitializer
        (
          "postprocessors",
          "algorithms applied to the argument after the optimization is concluded",
          io::Serialization::GetAgent( &m_PostProcessors),
          "()"
        );

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief performs pre-processing of the argument before the optimization is performed
      //! @param ARGUMENT data to be pre-processed
      //! @return pre-processed data
      virtual void PreProcess( t_ArgumentType &ARGUMENT) const
      {
      }

      //! @brief initiates the optimization
      //! @param ARGUMENT data to be optimized
      virtual void Optimize( t_ArgumentType &ARGUMENT) const = 0;

      //! @brief performs post-processing of the argument before the optimization is performed
      //! @param ARGUMENT data to be post-processed
      //! @return post-processed data
      virtual void PostProcess( t_ArgumentType &ARGUMENT) const
      {
      }

    ///////////////
    // operators //
    ///////////////

    public:

      //! @brief optimizes a given argument
      //! @param ARGUMENT data to be optimized
      void operator()( t_ArgumentType &ARGUMENT) const
      {
        // process and optimize the data
        PreProcess( ARGUMENT);
        for( auto proc_it( m_PreProcessors.Begin()), proc_it_end( m_PreProcessors.End()); proc_it != proc_it_end; ++proc_it)
        {
          ( **proc_it)( ARGUMENT);
        }
        Optimize( ARGUMENT);
        PostProcess( ARGUMENT);
        for( auto proc_it( m_PostProcessors.Begin()), proc_it_end( m_PostProcessors.End()); proc_it != proc_it_end; ++proc_it)
        {
          ( **proc_it)( ARGUMENT);
        }
      }

    }; // class OptimizationInterface

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_OPTIMIZATION_INTERFACE_H_
