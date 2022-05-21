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

#ifndef BCL_BIOL_BLAST_PROFILE_HANDLER_H_
#define BCL_BIOL_BLAST_PROFILE_HANDLER_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BlastProfileHandler
    //! @brief Handler class for BlastProfile object.
    //! @details Manages reading and writing BlastProfile objects from/to ProteinModel and AASequence
    //!
    //! @see @link example_biol_blast_profile_handler.cpp @endlink
    //! @author karakam
    //! @date 06/06/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BlastProfileHandler :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BlastProfileHandler();

      //! @brief virtual copy constructor
      BlastProfileHandler *Clone() const;

      //! @brief virtual destructor
      ~BlastProfileHandler();

    ///////////////////////
    // data access - Get //
    ///////////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read BLAST profile from ISTREAM for given amino acid
      //! @param ISTREAM input stream
      //! @param AMINO_ACID AABase into which blast profile is going to read
      //! @return std::istream which was read from
      static std::istream &ReadProfileForAA( std::istream &ISTREAM, AABase &AMINO_ACID);

      //! @brief write BLAST profile to OSTREAM for given amino acid
      //! @param OSTREAM output stream
      //! @param AMINO_ACID AABase from which blast profile is going to read
      //! @return std::ostream which was written to
      static std::ostream &WriteProfileForAA( std::ostream &OSTREAM, const AABase &AMINO_ACID);

      //! @brief read BLAST profile from ISTREAM for given AASequence
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which blast profile is going to read
      //! @return std::istream which was read from
      static std::istream &ReadProfileForAASequence( std::istream &ISTREAM, AASequence &AA_SEQUENCE);

      //! @brief write BLAST profile to OSTREAM for given AASequence
      //! @param OSTREAM output stream
      //! @param AA_SEQUENCE AASequence from which blast profile is going to read
      //! @return std::ostream which was written to
      static std::ostream &WriteProfileForAASequence( std::ostream &OSTREAM, const AASequence &AA_SEQUENCE);

      //! @brief read BLAST profile from ISTREAM for given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel into which blast profile is going to read
      //! @param PREFIX prefix of the protein, including path ( usually pdb id)
      //! @param EXTENSION the desired extension; if blank, .ascii6 will be preferred, but .ascii will also be accepted
      //! @return whether reading was successful
      static bool TryReadProfileForProteinModel
      (
        assemble::ProteinModel &PROTEIN_MODEL,
        const std::string &PREFIX,
        const std::string &EXTENSION = ""
      );

      //! @brief read all the blast profiles for a given sequence into a series of blast profiles
      //!        This is the only method of reading that does not change the underlying data structures and is likewise
      //!        safe inside threaded code that may be operating on the same sequence
      //! @param SEQUENCE the sequence of interest; only used to check that the blast profile has the same sequence
      //! @param PREFIX prefix of the protein, including path ( usually pdb id)
      //! @param EXTENSION the desired extension; if blank, .ascii6 will be preferred, but .ascii will also be accepted
      static storage::Vector< BlastProfile> ReadProfilesForConstAASequence
      (
        const AASequence &SEQUENCE,
        const std::string &PREFIX,
        const std::string &EXTENSION = ""
      );

      //! @brief read all the blast profiles for a given sequence into a series of blast profiles
      //!        This is the only method of reading that does not change the underlying data structures and is likewise
      //!        safe inside threaded code that may be operating on the same sequence
      //! @param SEQUENCE the sequence of interest; only used to check that the blast profile has the same sequence
      //! @param PROFILES reference to a vector that will be used to store the profiles
      //! @param ISTREAM input stream for reading the blast profiles
      static void ReadProfilesForConstAASequence
      (
        const AASequence &SEQUENCE,
        storage::Vector< BlastProfile> &PROFILES,
        std::istream &ISTREAM
      );

      //! @brief write BLAST profile to OSTREAM for given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel from which blast profile is going to read
      //! @param PREFIX prefix of the protein ( usually pdb id)
      //! @param PATH the directory path where blast files will be written
      //! @return whether writing was successful
      static bool WriteProfileForProteinModel
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const std::string &PREFIX,
        const std::string &PATH = "."
      );

    protected:

      //! @brief reads blast
      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; //class BlastProfile

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_BLAST_PROFILE_HANDLER_H_
