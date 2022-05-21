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

#ifndef BCL_CHEMISTRY_MOLECULE_ENSEMBLE_H_
#define BCL_CHEMISTRY_MOLECULE_ENSEMBLE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_molecule_complete.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeEnsemble
    //! @brief Container class for SmallMolecule objects
    //!
    //! @see @link example_chemistry_molecule_ensemble.cpp @endlink
    //! @author mendenjl
    //! @date Mar 08, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeEnsemble :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! list of molecules
      storage::List< MoleculeComplete> m_MoleculeEnsemble;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeEnsemble()
      {}

      //! @brief construct from list of SmallMoleculeConformations
      MoleculeEnsemble( const storage::List< MoleculeComplete> &MOLECULES) :
        m_MoleculeEnsemble( MOLECULES)
      {}

      //! @brief construct from an input stream linked to SDF format
      //! @param ISTREAM input stream, reads in SDF format
      //! @param RANGE a range of small molecules to load, by default loads all small molecules in the stream
      explicit MoleculeEnsemble
      (
        std::istream &ISTREAM,
        const math::Range< size_t> &RANGE = math::Range< size_t>( 0, std::numeric_limits< size_t>::max())
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeEnsemble
      MoleculeEnsemble *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return a shared pointer vector containing the ensemble
      const storage::List< MoleculeComplete> &GetMolecules() const
      {
        return m_MoleculeEnsemble;
      }

      //! @brief return a shared pointer vector containing the ensemble
      storage::List< MoleculeComplete> &GetMolecules()
      {
        return m_MoleculeEnsemble;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add a MoleculeComplete to the ensemble
      //! @param MOL_OBJ the object to add to the ensemble
      void PushBack( const MoleculeComplete &MOL_OBJ)
      {
        m_MoleculeEnsemble.PushBack( MOL_OBJ);
      }

      //! @brief merge this ensemble with another ensemble
      //! @param ENSEMBLE the object to add to the ensemble
      void Append( const MoleculeEnsemble &ENSEMBLE)
      {
        m_MoleculeEnsemble.Append( ENSEMBLE.m_MoleculeEnsemble);
      }

      //! @brief read additional molecules into the ensemble from an mdl stream
      //! @param ISTREAM input stream, reads in SDF format
      //! @param RANGE a range of small molecules to load, by default loads all small molecules in the stream
      void ReadMoreFromMdl
      (
        std::istream &ISTREAM,
        const math::Range< size_t> &RANGE = math::Range< size_t>( 0, std::numeric_limits< size_t>::max())
      );

      //! @brief GetSize returns size of the container
      size_t GetSize() const
      {
        return m_MoleculeEnsemble.GetSize();
      }

      //! @brief const iterator begin
      storage::List< MoleculeComplete>::const_iterator Begin() const
      {
        return m_MoleculeEnsemble.Begin();
      }
      //! @brief const iterator end
      storage::List< MoleculeComplete>::const_iterator End() const
      {
        return m_MoleculeEnsemble.End();
      }
      //! @brief iterator begin
      storage::List< MoleculeComplete>::iterator Begin()
      {
        return m_MoleculeEnsemble.Begin();
      }
      //! @brief iterator end
      storage::List< MoleculeComplete>::iterator End()
      {
        return m_MoleculeEnsemble.End();
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &WriteMDL( std::ostream &OSTREAM) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class MoleculeEnsemble

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MOLECULE_ENSEMBLE_H_
