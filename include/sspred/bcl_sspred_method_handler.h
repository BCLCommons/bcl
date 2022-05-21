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

#ifndef BCL_SSPRED_METHOD_HANDLER_H_
#define BCL_SSPRED_METHOD_HANDLER_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "io/bcl_io.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_environment_types.h"
#include "biol/bcl_biol_ss_types.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MethodHandler
    //! @brief the interface class for all secondary structure methods to be derived from.
    //! @details This class provided the interface for all various secondary structure prediction method classes to be derived
    //! from. It has no data storage, and it provides the interface for functions for getting three state or nine state
    //! predictions and reading predictions for a single amino acid and for a full sequence. These functions should be
    //! overwritten by the derived methods classes since each secondary structure method has their own format which
    //! requires distinct extracting, reading and writing
    //! functions.
    //!
    //! @see @link example_sspred_method_handler.cpp @endlink
    //! @author karakam
    //! @date Jun 3, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MethodHandler :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////
    // data //
    //////////

       //! enumerator for ssprediction output style
       enum OutputTypes
       {
         e_ThreeState,
         e_NineState,
         e_OriginalState
       };

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default cosntructor
      MethodHandler();

      //! @brief Clone function
      //! @return pointer to new MethodHandler
      MethodHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

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
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AMINO_ACID AABase into which secondary structure predictions will be read
      //! @param SS_METHOD Method to be read
      //! @return std::istream which was read from
      static std::istream &ReadPredictionsForAA
      (
        std::istream &ISTREAM,
        biol::AABase &AMINO_ACID,
        const Method &SS_METHOD
      );

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which secondary structure predictions will be read
      //! @param SS_METHOD Method to be read
      //! @return std::istream which was read from
      static std::istream &ReadPredictionsForAASequence
      (
        std::istream &ISTREAM,
        biol::AASequence &AA_SEQUENCE,
        const Method &SS_METHOD
      );

      //! @brief create all possible directory entries, a file could be read from and that exist
      //! @param SS_METHOD Method to be read
      //! @param CHAIN_ID possible chainid to add to prefix, so that file can be found
      //! @param PREFIX prefix of the sequence ( usually pdb id)
      //! @param PATH path where secondary structure prediction files can be found
      //! @return vector of directory entries, that actually exist
      static storage::Vector< io::DirectoryEntry> PossibleFileEntries
      (
        const Method &SS_METHOD,
        const char CHAIN_ID,
        const std::string &PREFIX,
        const std::string &PATH = "."
      );

      //! @brief create all possible directory entries, a file could be read from and that exist
      //! @param SS_METHOD Method to be read
      //! @param EXTENSION extension of the ss pred file
      //! @param CHAIN_ID possible chainid to add to prefix, so that file can be found
      //! @param PREFIX prefix of the sequence ( usually pdb id)
      //! @param PATH path where secondary structure prediction files can be found
      //! @return vector of directory entries, that actually exist
      static storage::Vector< io::DirectoryEntry> PossibleFileEntries
      (
        const std::string &EXTENSION,
        const char CHAIN_ID,
        const std::string &PREFIX,
        const std::string &PATH
      );

      //! @brief create a subset of the methods given, for which the actual files exist
      //! If PREFIX of "1UBI" is supplied , the function will automatically supply the ChainID if necessary to create "1UBI_"
      //! and with the SSMethod, the right extension will be supplied leading to "1UBI_.jufo", supplying an input
      //! path will help if the file is somewhere else than the the executable like: "/home/user/" leading to
      //! "/home/user/1UBI_.jufo" - if any of the files exist, the method will be added to the set returned
      //! @param SS_METHODS Set of Method to be read
      //! @param CHAIN_ID possible chainid to add to prefix, so that file can be found
      //! @param PREFIX prefix of the sequence ( usually pdb id)
      //! @param PATH path where secondary structure prediction files can be found
      //! @return subset of methods, for which files exist
      static storage::Set< Method> AvailablePredictionFiles
      (
        const storage::Set< Method> &SS_METHODS,
        const char CHAIN_ID,
        const std::string &PREFIX,
        const std::string &PATH = "."
      );

      //! @brief read specified set of METHODS for AA_SEQUENCE
      //! If PREFIX of "1UBI" is supplied , the function will automatically supply the ChainID leading to "1UBI_"
      //! and with the SSMethod, the right extension will be supplied leading to "1UBI_.jufo", supplying an input
      //! path will help if the file is somewhere else than the the executable like: "/home/user/" leading to
      //! "/home/user/1UBI_.jufo"
      //! @param SS_METHODS Set of Method to be read
      //! @param AA_SEQUENCE AASequence into which secondary structure predictions will be read
      //! @param PREFIX prefix of the sequence ( usually pdb id)
      //! @param PATH path where secondary structure prediction files can be found
      //! @return whether reading has been successful
      static bool ReadPredictionsForAASequence
      (
        const storage::Set< Method> &SS_METHODS,
        biol::AASequence &AA_SEQUENCE,
        const std::string &PREFIX,
        const std::string &PATH = "."
      );

      //! @brief read specified set of METHODS for PROTEIN_MODEL
      //! If PREFIX of "1UBI" is supplied , the function will automatically supply the ChainID leading to "1UBI_"
      //! and with the SSMethod, the right extension will be supplied leading to "1UBI_.jufo", supplying an input
      //! path will help if the file is somewhere else than the the executable like: "/home/user/" leading to
      //! "/home/user/1UBI_.jufo"
      //! @param SS_METHODS Set of Method to be read
      //! @param PROTEIN_MODEL ProteinModel into which secondary structure predictions will be read
      //! @param PREFIX prefix of the sequence ( usually pdb id)
      //! @param PATH path where secondary structure prediction files can be found
      //! @return whether reading has been successful
      static bool ReadPredictionsForProteinModel
      (
        const storage::Set< Method> &SS_METHODS,
        assemble::ProteinModel &PROTEIN_MODEL,
        const std::string &PREFIX,
        const std::string &PATH = "."
      );

      //! @brief read all SSPred methods for which prediction files are available for a given protein model
      //! If PREFIX of "1UBI" is supplied , the function will automatically supply the ChainID leading to "1UBI_"
      //! and with the SSMethod, the right extension will be supplied leading to "1UBI_.jufo", supplying an input
      //! path will help if the file is somewhere else than the the executable like: "/home/user/" leading to
      //! "/home/user/1UBI_.jufo"
      //! @param PROTEIN_MODEL ProteinModel into which secondary structure predictions will be read
      //! @param PREFIX prefix of the sequence ( usually pdb id)
      //! @param PATH path where secondary structure prediction files can be found
      static void ReadAllPredictionsForProteinModel
      (
        assemble::ProteinModel &PROTEIN_MODEL,
        const std::string &PREFIX,
        const std::string &PATH = "."
      );

      //! @brief write three state secondary structure predictions for given sequence from the provided OSTREAM
      //! @param OSTREAM output stream
      //! @param AA_SEQUENCE AASequence for which secondary structure predictions will be written
      //! @param SS_METHOD Method to be written
      //! @param OUTPUT_TYPE enum to decide the output_type, set to three-state by default
      //! @return std::ostream which was written to
      static std::ostream &WritePredictionsForAASequence
      (
        std::ostream &OSTREAM,
        const biol::AASequence &AA_SEQUENCE,
        const Method &SS_METHOD,
        const OutputTypes OUTPUT_TYPE = e_ThreeState
      );

      //! @brief writes out the specified SS_METHODS to files with given PREFIX to provided path for given AA_SEQUENCE
      //! this function checks that given SS_METHODS are defined for this sequence and secondly that
      //! they are not overwritten in the file system. output is directed to std::cout if WRITE_TO_STD_COUT is set
      //! @param SS_METHODS vector of SSMethods to be read
      //! @param AA_SEQUENCE AASequence for which secondary structure predictions will be written
      //! @param PREFIX prefix of the SSPredictions files
      //! @param WRITE_TO_STD_COUT boolean to whether direct output to std::cout
      //! @param PATH directory path where the files for SSPredictions reside
      //! @param OUTPUT_TYPE enum to decide the output_type, set to three-state by default
      //! @return whether writing was successful
      static bool WritePredictionsForAASequence
      (
        const storage::Set< Method> &SS_METHODS,
        const biol::AASequence &AA_SEQUENCE,
        const std::string &PREFIX,
        const bool &WRITE_TO_STD_COUT,
        const std::string &PATH = ".",
        const OutputTypes OUTPUT_TYPE = e_ThreeState
      );

      //! @brief writes out the specified SS_METHODS to files with given PREFIX to provided path for given PROTEIN_MODEL
      //! this function checks that given SS_METHODS are defined for this sequence and secondly that
      //! they are not overwritten in the file system. output is directed to std::cout if WRITE_TO_STD_COUT is set
      //! @param SS_METHODS vector of SSMethods to be read
      //! @param PROTEIN_MODEL ProteinModel for which secondary structure predictions will be written
      //! @param PREFIX prefix of the SSPredictions files
      //! @param WRITE_TO_STD_COUT boolean to whether direct output to std::cout
      //! @param PATH directory path where the files for SSPredictions reside
      //! @param OUTPUT_TYPE enum to decide the output_type, set to three-state by default
      //! @return whether writing was successful
      static bool WritePredictionsForProteinModel
      (
        const storage::Set< Method> &SS_METHODS,
        const assemble::ProteinModel &PROTEIN_MODEL,
        const std::string &PREFIX,
        const bool &WRITE_TO_STD_COUT,
        const std::string &PATH = ".",
        const OutputTypes OUTPUT_TYPE = e_ThreeState
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief initializes predictions of all amino acids in AA_SEQUENCE for SS_METHOD to provided PREDICTION
      //! @param SS_METHOD Method of interest
      //! @param AA_SEQUENCE AASequence to be updated
      //! @param PREDICTION MethodInterface to be assigned
      static void InitializePredictionsForAASequence
      (
        const Method &SS_METHOD,
        biol::AASequence &AA_SEQUENCE,
        const MethodInterface &PREDICTION
      );

      //! @brief function for predictions for given SS_METHOD in given subsequence intervals to provided PREDICTION
      //! @param SS_METHOD Method of interest
      //! @param AA_SEQUENCE AASequence to be updated
      //! @param PREDICTION MethodInterface to be assigned
      //! @param FIRST index of the first amino acid of subsequence
      //! @param LAST index of the last amino acid of of subsequence
      static void SetPredictionsForSubSequence
      (
        const Method &SS_METHOD,
        biol::AASequence &AA_SEQUENCE,
        const MethodInterface &PREDICTION,
        const size_t FIRST,
        const size_t LAST
      );

    private:

    }; // class MethodHandler

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_METHOD_HANDLER_H_
