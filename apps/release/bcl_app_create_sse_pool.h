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

#ifndef BCL_APP_CREATE_SSE_POOL_H_
#define BCL_APP_CREATE_SSE_POOL_H_

// include header of this class
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "storage/bcl_storage_map.h"

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CreateSSEPool
    //! @brief This application generates the SSEPool for the given proteins
    //! @details This class generates the secondary structure pool for given proteins using individual and consensus
    //! secondary structure predictions
    //!
    //! @see @link example_app_create_sse_pool.cpp @endlink
    //! @author karakam
    //! @date 03/24/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CreateSSEPool :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! fasta tag
      util::ShPtr< command::FlagStatic> m_PrefixFlag;

      //! output path flag
      util::ShPtr< command::FlagStatic> m_OutputPrefixFlag;

      //! chain id
      util::ShPtr< command::FlagStatic> m_ChainIdFlag;

      //! pdb tag
      util::ShPtr< command::FlagStatic> m_PdbFlag;
      util::ShPtr< command::Parameter> m_PdbPoolFileOutParam;

      //! evaluate_pool flag
      util::ShPtr< command::FlagStatic> m_EvaluatePoolFlag;

      //! ss methods to be used in pool generation
      util::ShPtr< command::FlagInterface> m_SsMethods;

      //! flag to indicate whether long sses should be chopped into two parts
      util::ShPtr< command::FlagInterface> m_ChopSsesFlag;

      //! flag to indicate whether user wants to add systematically extended copies of the sses in the pool
      util::ShPtr< command::FlagInterface> m_SystematicExtendFlag;

      //! threshold for defining sse regions
      util::ShPtr< command::FlagStatic> m_SseThresholdFlag;
      util::ShPtr< command::ParameterInterface> m_HelixThresholdParam;
      util::ShPtr< command::ParameterInterface> m_StrandThresholdParam;
      util::ShPtr< command::ParameterInterface> m_CoilThresholdParam;

      //! flag for factory
      util::ShPtr< command::FlagStatic> m_FactoryFlag;
      util::ShPtr< command::Parameter>  m_FactoryParam;
      util::ShPtr< command::Parameter>  m_FactoryStreamParam;

      //! flag for joining and separating adjacent sses of the same type
      util::ShPtr< command::FlagStatic> m_JoinSeparateFlag;

      //! thresholds map for prediction thresholds
      mutable storage::Map< biol::SSType, double> m_ThresholdsMap;

      //! min size map for sses in the pool
      mutable storage::Map< biol::SSType, size_t> m_SSEMinSizeMap;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      CreateSSEPool();

    public:

      static const ApplicationType CreateSSEPool_Instance;

      //! @brief Clone function
      //! @return pointer to new CreateSSEPool
      CreateSSEPool *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief process pool by joining, separating, chopping ...
      //! @param POOL the pool to process
      void ProcessPool( assemble::SSEPool &POOL) const;

    }; // class CreateSSEPool
  } // namespace app
} // namespace bcl
#endif // BCL_APP_CREATE_SSE_POOL_H_
