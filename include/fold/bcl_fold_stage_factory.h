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

#ifndef BCL_FOLD_STAGE_FACTORY_H_
#define BCL_FOLD_STAGE_FACTORY_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_stage_interface.h"
#include "mc/bcl_mc_stage.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StageFactory
    //! @brief Factory class for producing stages for fold minimization.
    //! @details StageFactory class allows creating stages for fold minimization. If a single stage, default behavior,
    //! fold minimization is requested then it will construct it using the command line flags. Otherwise, it will use
    //! the command line provided stage file to construct the stages. This class is responsible for collecting all
    //! different protocols used in all stages and calling the functions for initializing scores and mutates.
    //!
    //! @remarks example unnecessary
    //! @author karakam, fischea
    //! @date Feb 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StageFactory :
      public util::ObjectInterface
    {
    public:

    //////////
    // data //
    //////////

      //! enumerator for line types
      enum LineType
      {
        e_NumberCycles,         // Number of cycles
        e_StageStart,           // Start of a stage definition, followed by an optional name
        e_Type,                 // type of the approximation used for this stage
        e_FoldProtocols,        // fold protocols to use
        e_ScoreProtocols,       // score protocols to use (optional)
        e_MutateProtocols,      // mutate protocols to use (optional)
        e_ScoreWeightSet,       // score weight set table
        e_ScoreWeightSetFile,   // file containing score weight set
        e_ScoreDropoutRate,     // Fraction of scores to drop randomly for each protein model
        e_MutateWeightSet,      // mutate weight set table
        e_MutateWeightSetFile,  // file containing mutate weighStaget set
        e_NumberIterations,     // followed by total number iterations and max number unimproved
        e_ModifyStartModel,     // true if the stage should modify the start model, false otherwise
        e_PrintStartModel,      // whether to print the start model for this stage
        e_PrintIterationModel,  // whether to print the every iteration model for this stage
        e_PrintEndModel,        // whether to print the end model for this stage
        e_PrintTrackerHistory,  // whether to print the tracker history for this stage
        e_PoolPostfix,          // file postfix containing pool that should be used for this stage
        e_StageEnd,             // end of a stage definition
        s_NumberLineTypes
      };

      //! line type names
      static const std::string s_LineTypeNames[];

      //! @brief finds the LineType enum that corresponds to given string
      //! @param LINE_NAME Line name of interest
      //! @return the LineType enum that corresponds to given string
      static LineType LineTypeFromString( const std::string &LINE_NAME);

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      StageFactory();

      //! @brief Clone function
      //! @return pointer to new StageFactory
      StageFactory *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      ////////////////stages_cycles
      // operations //
      ////////////////

      //! @brief construct stages
      //! @return vector of stages
      static util::ShPtrVector< StageInterface> CreateStages();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief construct and return stages from stage file
      //! @param NUMBER_CYCLES number of cycles which can be updated from command line
      //! @return stages constructed from stage file
      static util::ShPtrVector< StageInterface> CreateStagesFromFile( size_t &NUMBER_CYCLES);

      //! @brief constructs one stage from an IFStream
      //! @param READ IFStream to create the stage from
      //! @param NAME the name of the stage
      //! @param NUMBER the number of the stage
      //! @return shared pointer to the constructed stage
      static util::ShPtr< StageInterface> CreateStage( io::IFStream &READ, const std::string &NAME, const size_t &NUMBER);

      //! @brief constructs a stage from the command line
      //! @return shared pointer to the constructed stage
      static util::ShPtr< StageInterface> CreateStage();

      //! @brief creates a monte carlo metropolis stage based on the given parameters
      //! @param SETTINGS contains the parameters to construct the stage from
      //! @return shared pointer to the constructed stage
      static util::ShPtr< StageInterface> CreateMcStage
      (
        storage::HashMap< size_t, storage::Vector< std::string> > &SETTINGS,
        const size_t &NUMBER
      );

      //! @brief sets the temperature for the given STAGE
      static void SetTemperature( mc::Stage &STAGE);

      //! @brief write stage to stream
      //! @param STAGE Stage to be written
      //! @param OSTREAM ostream to be written to
      //! @return ostream which was written to
      static std::ostream &WriteStage( const mc::Stage &STAGE, std::ostream &OSTREAM);

      //! @brief write given protocol list to stream
      //! @param PROTOCOL_LIST Protocol list of interest
      //! @param OSTREAM ostream to be written to
      //! @return ostream which was written to
      static std::ostream &WriteProtocolList( const storage::List< Protocol> &PROTOCOL_LIST, std::ostream &OSTREAM);

      //! @brief removes duplicates if any from the given protocol list
      //! @param PROTOCOL_LIST Protocol list of interest
      //! @return list of unique protocols
      static storage::List< Protocol> CreateUniqueProtocolList( const storage::List< Protocol> &PROTOCOL_LIST);

      //! @brief construct a protocol list from the given vector of strings making sure there are no duplicates
      //! @param STRING_VECTOR Vector of strings
      //! @return list of unique protocols
      static storage::List< Protocol> ConstructProtocolList( const storage::Vector< std::string> &STRING_VECTOR);

      //! @brief read table from file
      //! @param FILENAME name for the file that contains table
      //! @return Table read from the file
      static storage::Table< double> ReadTableFromFile( const std::string &FILENAME);

      //! @brief write table to file
      //! @param TABLE Table to write
      //! @param FILENAME name for the file to which the table will be written to
      static void WriteTableToFile( const storage::Table< double> &TABLE, const std::string &FILENAME);

      //! @brief function to set the scoring function for the given single stage
      //! @param SINGLE_STAGE Stage of interest
      static void SetSingleStageScoreFunction( mc::Stage &SINGLE_STAGE);

      //! @brief function to set the scoring weight set for the given single stage
      //! @param SINGLE_STAGE Stage of interest
      static void SetSingleStageScoreWeightSet( mc::Stage &SINGLE_STAGE);

      //! @brief function to set the mutate for the given single stage
      //! @param SINGLE_STAGE Stage of interest
      static void SetSingleStageMutate( mc::Stage &SINGLE_STAGE);

      //! @brief function to set the mutate tree for the given single stage
      //! @param SINGLE_STAGE Stage of interest
      static void SetSingleStageMutateTree( mc::Stage &SINGLE_STAGE);

    }; // class StageFactory

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_STAGE_FACTORY_H_
