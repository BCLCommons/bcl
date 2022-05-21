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

#ifndef BCL_PDB_LIGAND_H_
#define BCL_PDB_LIGAND_H_

// include the namespace header
#include "bcl_pdb.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_pdb_line.h"
#include "bcl_pdb_residue.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Ligand
    //! @brief Ligand is a pdb residue with non standard residue type
    //! @details besides the normal residue information, it has the full name, formula and the atom connectivities
    //!
    //! @see @link example_pdb_ligand.cpp @endlink
    //! @author woetzen, alexanns, loweew
    //! @date Jun 24, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Ligand :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! full name from HETNAM
      std::string m_Fullname;

      //! formula from FORMUL lines
      std::string m_Formula;

      std::string            m_ResidueName; //! residue name
      char                   m_ChainID;     //! chain ID
      int                    m_PDBID;       //! pdb id
      char                   m_ICode;       //! insertion code
      util::ShPtrList< Line> m_Lines;       //! pdb lines as sharedpointervector

      //! the connections within the ligand
      storage::Map< size_t, storage::Set< size_t> > m_InternalConnections;

      //! the connections outside the ligand
      storage::Map< size_t, storage::Set< size_t> > m_ExternalConnections;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Ligand();

      //! @brief constructor from lines
      Ligand( const util::ShPtrList< Line> &LINES);

      //! @brief Clone function
      //! @return pointer to new Ligand
      Ligand *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! return ID
      int const &GetPDBID() const
      {
        return m_PDBID;
      }

      //! set ID
      void SetPDBID( const int PDBID)
      {
        m_PDBID = PDBID;
      }

      //! return insertion code
      char const &GetICode() const
      {
        return m_ICode;
      }

      //! set insertion code
      void SetICode( const char ICODE)
      {
        m_ICode = ICODE;
      }

      //! return chain id
      char const &GetChainID() const
      {
        return m_ChainID;
      }

      //! set chain id
      void SetChainID( const char CHAINID)
      {
        m_ChainID = CHAINID;
      }

      //! return residue name
      const std::string &GetResidueName() const
      {
        return m_ResidueName;
      }

      //! set residue name
      void SetResidueName( const std::string &RESIDUENAME)
      {
        m_ResidueName = RESIDUENAME;
      }

      //! return Lines
      util::ShPtrList< Line> const &GetLines() const;

      //! change Lines
      util::ShPtrList< Line> &ChangeLines();

      //! @brief access to the fullname
      //! @return the full ligand name
      const std::string &GetFullname() const;

      //! @brief access to the formula
      //! @return the ligand formula
      const std::string &GetFormula() const;

      //! @brief set full name
      //! @param FULL_NAME as reported for that residue name in the HETNAME line
      void SetFullname( const std::string &FULL_NAME);

      //! @brief set the formula
      //! @param FORMULA as reported in the FORMUL lines for the residue name
      void SetFormula( const std::string &FORMULA);

      //! @brief add connections for an atom serial
      //! @param CONNECTIONS map of atom serials with sets of all atoms that are connected; all connections are present
      //!        twice and inverted, as they come from the pdb
      void AddConnections( const storage::Map< size_t, storage::Set< size_t> > &CONNECTIONS);

    ////////////////
    // operations //
    ////////////////

      //! @brief get all atom serials for that residue
      //! @return set of atom serials
      storage::Set< size_t> AtomSerials() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class Ligand

  } // namespace pdb
} // namespace bcl

#endif // BCL_PDB_LIGAND_H_ 
