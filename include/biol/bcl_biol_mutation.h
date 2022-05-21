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

#ifndef BCL_BIOL_MUTATION_H_
#define BCL_BIOL_MUTATION_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_base.h"
#include "bcl_biol_aa_types.h"
#include "io/bcl_io_serialization_interface.h"
#include "util/bcl_util_object_instances.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Mutation
    //! @brief Mutation of a particular protein or sequence
    //! @see @link example_biol_mutation.cpp @endlink
    //! @author mendenjl
    //! @date Jan 17, 2019
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Mutation :
      public io::SerializationInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @brief Residue ID
      int m_ResidueNumber;

      //! @brief native residue type
      AAType m_NativeType;

      //! @brief Mutant residue type
      AAType m_MutantType;

      //! @brief simple pointer to the mutated residue(s)
      mutable util::SiPtrVector< const AABase> m_AAs;

      //! chain ids; blank if chain ids aren't known yet
      mutable std::string m_ChainIDs;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Mutation() :
        m_ResidueNumber( 0)
      {
      }

      //! @brief constructor from given input data
      //! @param RESNUM - The residue number in the PDB
      //! @param NATIVE - The native residue type
      //! @param MUTANT - The mutant residue type
      Mutation( const int RESNUM, const AAType &NATIVE, const AAType &MUTANT, const std::string &CHAIN_IDS = std::string()) :
        m_ResidueNumber( RESNUM),
        m_NativeType( NATIVE),
        m_MutantType( MUTANT),
        m_ChainIDs( CHAIN_IDS)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Mutation
      Mutation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief Get residue number
      //! @return the residue number
      const int &GetResidueNumber() const
      {
        return m_ResidueNumber;
      }

      //! @brief Get the native type
      //! @return The native type
      const AAType &GetNativeType() const
      {
        return m_NativeType;
      }

      //! @brief Get the mutant type
      //! @return the mutant type
      const AAType &GetMutantType() const
      {
        return m_MutantType;
      }

      //! @brief Get the chain ids
      //! @return the chain ids
      const std::string &GetChainIDs() const
      {
        return m_ChainIDs;
      }

      //! @brief write the mutation as a string in standard mutation format
      std::string ToString() const;

      //! @brief write the mutation as a string in standard mutation format
      static Mutation FromString( const std::string &STR);

      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      util::ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const;

      //! @brief connect the mutation with a particular residue
      const util::SiPtrVector< const AABase> &GetAAs() const
      {
        return m_AAs;
      }

      //! @brief connect the mutation with a particular residue
      void SetAA( const AABase &BASE) const;

    ////////////////
    // operators  //
    ////////////////

      //! @brief compare two Mutation
      //! @param POINT the point to compare to this point
      //! @return if *this and POINT have identical data
      bool operator ==( const Mutation &POINT) const
      {
        return
            m_ResidueNumber == POINT.m_ResidueNumber &&
            m_NativeType == POINT.m_NativeType &&
            m_MutantType == POINT.m_MutantType &&
            m_ChainIDs == POINT.m_ChainIDs;
      }

      //! @brief compare two Mutation
      //! @param POINT the point to compare to this point
      //! @return if *this and POINT have identical data
      bool operator <( const Mutation &POINT) const
      {
        return
            m_ResidueNumber < POINT.m_ResidueNumber
            ||
            (
              m_ResidueNumber == POINT.m_ResidueNumber
              && m_NativeType == POINT.m_NativeType
              && m_MutantType < POINT.m_MutantType
            );
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryRead( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

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

    }; // class Mutation
  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_MUTATION_H_
