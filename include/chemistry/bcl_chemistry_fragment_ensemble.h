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

#ifndef BCL_CHEMISTRY_FRAGMENT_ENSEMBLE_H_
#define BCL_CHEMISTRY_FRAGMENT_ENSEMBLE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_hydrogens_handler.h"
#include "math/bcl_math_range.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentEnsemble
    //! @brief Container class for SmallMolecule objects
    //!
    //! @see @link example_chemistry_fragment_ensemble.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Feb 21, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentEnsemble :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! list of share pointers to SmallMoleculeConformations
      storage::List< FragmentComplete> m_FragmentEnsemble;

    public:

    //////////////
    // typedefs //
    //////////////

      typedef storage::List< FragmentComplete>::const_iterator const_iterator;
      typedef storage::List< FragmentComplete>::iterator iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentEnsemble()
      {}

      //! @brief construct from list of SmallMoleculeConformations
      FragmentEnsemble( const storage::List< FragmentComplete> &MOLECULES) :
        m_FragmentEnsemble( MOLECULES)
      {}

      //! @brief construct from an input stream linked to SDF format
      //! @param ISTREAM input stream, reads in SDF format
      //! @param ADD_H whether to add hydrogens to the ensemble automatically if the format is e_BCL
      //! @param RANGE a range of small molecules to load, by default loads all small molecules in the stream
      explicit FragmentEnsemble
      (
        std::istream &ISTREAM,
        const sdf::HydrogenHandlingPref &H_PREF = sdf::e_Maintain,
        const math::Range< size_t> &RANGE = math::Range< size_t>( 0, std::numeric_limits< size_t>::max()),
        const sdf::NeutralizationPref &NEUTRALIZATION = sdf::e_CmdLine
      );

      //! @brief Clone function
      //! @return pointer to new FragmentEnsemble
      FragmentEnsemble *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return a shared pointer vector containing the ensemble
      const storage::List< FragmentComplete> &GetMolecules() const
      {
        return m_FragmentEnsemble;
      }

      //! @brief return a shared pointer vector containing the ensemble
      storage::List< FragmentComplete> &GetMolecules()
      {
        return m_FragmentEnsemble;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add a FragmentComplete to the ensemble
      //! @param MOL_OBJ the object to add to the ensemble
      void PushBack( const FragmentComplete &MOL_OBJ)
      {
        m_FragmentEnsemble.PushBack( MOL_OBJ);
      }

      //! @brief remove a FragmentComplete to the ensemble
      void PopBack()
      {
        m_FragmentEnsemble.PopBack();
      }

      //! @brief merge this ensemble with another ensemble
      //! @param ENSEMBLE the object to add to the ensemble
      void Append( const FragmentEnsemble &ENSEMBLE)
      {
        m_FragmentEnsemble.Append( ENSEMBLE.m_FragmentEnsemble);
      }

      //! @brief shuffle the ensemble
      void Shuffle();

      //! @brief read additional molecules into the ensemble
      //! @param ISTREAM input stream, reads in SDF format
      //! @param ADD_H whether to add hydrogens to the ensemble automatically if the format is e_BCL
      //! @param RANGE a range of small molecules to load, by default loads all small molecules in the stream
      void ReadMoreFromMdl
      (
        std::istream &ISTREAM,
        const sdf::HydrogenHandlingPref &H_PREF = sdf::e_Maintain,
        const math::Range< size_t> &RANGE = math::Range< size_t>( 0, std::numeric_limits< size_t>::max()),
        const sdf::NeutralizationPref &NEUTRALIZATION = sdf::e_CmdLine
      );

      //! @brief GetSize returns size of the container
      size_t GetSize() const
      {
        return m_FragmentEnsemble.GetSize();
      }

      //! @brief return true if the ensemble is empty
      bool IsEmpty() const
      {
        return m_FragmentEnsemble.IsEmpty();
      }

      //! @brief const iterator begin
      const_iterator Begin() const
      {
        return m_FragmentEnsemble.Begin();
      }

      //! @brief const iterator end
      const_iterator End() const
      {
        return m_FragmentEnsemble.End();
      }

      //! @brief iterator begin
      iterator Begin()
      {
        return m_FragmentEnsemble.Begin();
      }

      //! @brief iterator end
      iterator End()
      {
        return m_FragmentEnsemble.End();
      }

      //! @brief Remove all hydrogens
      void RemoveH();

      //! @brief Add hydrogens
      void SaturateWithH();

      //! @brief sort ensemble based on increasing value property
      //! @param PROPERTY property of interest which is stored on the molecule
      //! @return a vector with molecules sorted on the basis of property in increasing order
      void Sort( const std::string PROPERTY_NAME);

      //! @brief sort ensemble based on increasing value property
      //! @param PROPERTY property of interest which is stored on the molecule
      //! @return a vector with molecules sorted on the basis of property in increasing order
      void Sort( const storage::Vector< std::string> PROPERTY_NAMES);

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

    }; // class FragmentEnsemble

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_ENSEMBLE_H_
