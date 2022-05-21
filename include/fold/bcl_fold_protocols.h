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

#ifndef BCL_FOLD_PROTOCOLS_H_
#define BCL_FOLD_PROTOCOLS_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_protocol_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Protocols
    //! @brief enumerator class for all different fold protocols
    //! @details this util::Enumerate derived class allows enumeration of various fold protocols. By default, it only
    //! defines the e_Default which is the default de novo folding protocol. Each protocol is responsible for adding
    //! itself to this enum. This way, only the protocols compiled are accessible through fold application
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Nov 11, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Protocols :
      public util::Enumerate< util::SiPtr< ProtocolInterface>, Protocols>
    {
      friend class util::Enumerate< util::SiPtr< ProtocolInterface>, Protocols>;

    public:

    //////////
    // data //
    //////////

      Protocol e_Default;                       //!< default protocol
      Protocol e_Create;                        //!< creation protocol
      Protocol e_Assembly;                      //!< assembly protocol
      Protocol e_Refinement;                    //!< refinement protocol
      Protocol e_Membrane;                      //!< membrane protocol
      Protocol e_Restraint;                     //!< restraint protocol
      Protocol e_EM;                            //!< EM protocol
      Protocol e_Multimer;                      //!< multimer protocol
      Protocol e_Template   ;                   //!< template protocol which assembles protein models from templates
      Protocol e_LoopCoordinateAdd;             //!< loop coordinate add protocol
      Protocol e_LoopClose;                     //!< loop close protocol
      Protocol e_Dock;                          //!< dock protocol
      Protocol e_Ensemble;                      //!< ensemble protocol
      Protocol e_EnsembleSwitchConformation;    //!< ensemble protocol for switching conformation
      Protocol e_EnsembleReplicateConformation; //!< ensemble protocol for replicating conformation
      Protocol e_EnsembleFilter;                //!< ensemble protocol for filtering conformations

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Protocols();

    public:

      //! @brief return command line flag for choosing the protocols to be used to run folding
      //! @return command line flag for choosing the protocols to be used to run folding
      static util::ShPtr< command::FlagInterface> &GetFlagProtocols();

      //! @brief return command line flag for choosing the protocols to be used for mutates
      //! @return command line flag for choosing the protocols to be used for mutates
      static util::ShPtr< command::FlagInterface> &GetFlagMutateProtocols();

      //! @brief return command line flag for choosing the protocols to be used for scoring
      //! @return command line flag for choosing the protocols to be used for scoring
      static util::ShPtr< command::FlagInterface> &GetFlagScoreProtocols();

      //! @brief return list of protocols specified to be used specified in command line
      //! @return list of protocols to be used specified in command line
      static storage::List< Protocol> GetCommandLineProtocolList();

      //! @brief return list of protocols to be used for defining scores specified in command line
      //! @return list of protocols to be used for defining scores specified in command line
      static storage::List< Protocol> GetCommandLineScoreProtocolList();

      //! @brief return list of protocols to be used for defining mutates specified in command line
      //! @return list of protocols to be used for defining mutates specified in command line
      static storage::List< Protocol> GetCommandLineMutateProtocolList();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class Protocols

    //! @brief construct on access function for all Protocols
    //! @return reference to only instances of Protocols
    BCL_API
    Protocols &GetProtocols();

  } // namespace fold

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< SiPtr< fold::ProtocolInterface>, fold::Protocols>;

  } // namespace util
} // namespace bcl

#endif // BCL_FOLD_PROTOCOLS_H_
