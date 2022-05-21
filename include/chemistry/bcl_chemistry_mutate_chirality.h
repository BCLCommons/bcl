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

#ifndef BCL_CHEMISTRY_MUTATE_CHIRALITY_H_
#define BCL_CHEMISTRY_MUTATE_CHIRALITY_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "graph/bcl_graph_const_graph.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateChirality
    //! @brief This class is for randomly using different rotamers of a fragment to change conformation of a molecule of interest.
    //! @details
    //!
    //! @see @link example_chemistry_mutate_chirality.cpp @endlink
    //! @author kothiwsk
    //! @date Mar 31, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateChirality :
      public math::MutateInterface< FragmentComplete>
    {

    private:

    //////////
    // data //
    //////////

      graph::ConstGraph< size_t, size_t>              m_Graph;
      graph::ConstGraph< size_t, size_t>              m_GraphBondOrder;

      storage::Vector< size_t>                        m_ChiralCenters;

      storage::Vector< size_t>                        m_UnknownChiralCenters;

      bool                                            m_ChangeUnknownChiralOnly;

      storage::Vector< storage::Map< size_t, storage::Vector< size_t> > > m_BondsToDownstreamAtoms;

      storage::Vector< size_t>                        m_OnlyHasIsometry;

      //! cached bond info
      storage::Vector< sdf::BondInfo>                 m_BondInfo;

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
      MutateChirality()
      {
      }

      //! @brief constructor taking the member variable parameters
      //! @param FRAGMENT fragment which the mutate mutates
      //! @param UNKOWN_CHIRAL_ONLY boolean to specify whether to change unknown chiral only
      //! @param MOLECULE molecule of interest
      MutateChirality( const FragmentComplete &FRAGMENT, bool UNKNOWN_CHIRAL_ONLY = true);

      //! @brief Clone function
      //! @return pointer to new MutateChirality
      MutateChirality *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns atom indices that are connected to different bonds
      //! @return atom indices that are connected to different rotatable bonds
      const storage::Vector< size_t> &GetChiralCenters() const;

      //! @brief returns bond info (cached for performance)
      //! @return bond info
      const storage::Vector< sdf::BondInfo> &GetBondInfo() const
      {
        return m_BondInfo;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a conformation and returns a mutated conformation
      //! @param MOLECULE conformation of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< FragmentComplete> operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief rotate the molecule of interest at a particular bond at a particular angle
      //! @param ATOM_INFO the atom info that needs to be updated with coordinates after rotating about a particular bond
      //! @param CONNECTED_ATOMS the connected atoms associated with bond of interest
      //! @param NEW_ANGLE the new angle for the bond of interest
      //! @param EXISTING_ANGLE the angle for the bond of interest in the argument that is passed in
      //! @param FORCE_ALTER_LABEL alter the chirality label even if the chirality cannot be adjusted via reflections
      //!        This is useful if later code has to assume that the mutation was successful
      FragmentComplete MutateSpecificChirality
      (
        const FragmentComplete &MOLECULE,
        const size_t &ATOM_INDEX,
        const bool &ALTER_CHIRALITY_LABEL = true,
        const bool &FORCE_ALTER_LABEL = false
      ) const;

      //! @brief Change the chiral, pseudochiral, or isometric-bond-adjacent identity of a given atom
      //! @param MOLECULE the molecule of interest
      //! @param ATOM_INDEX the atom whose chirality or E/Z bond identity should be updated
      //! @param ALTER_CHIRALITY_LABEL whether to swap/update the chirality or isometry label
      //! @param FORCE_ALTER_LABEL alter the chirality label even if the chirality cannot be adjusted via reflections
      //!        This is useful if later code has to assume that the mutation was successful
      bool MutateSpecificChirality
      (
        AtomVector< AtomComplete> &MOLECULE,
        const size_t &ATOM_INDEX,
        const bool &ALTER_CHIRALITY_LABEL = true,
        const bool &FORCE_ALTER_LABEL = false
      ) const;

      //! @brief Copy the chirality of a different molecule onto this molecule
      //! This function assumes that MOLECULE currently has up-to-date chirality information
      //! Atoms with ambiguous or unknown chirality from TEMPLATE are ignored
      //! @param FORCE_ALTER_LABEL alter the chirality label even if the chirality cannot be adjusted via reflections
      //!        This is useful if later code has to assume that the mutation was successful
      //! @return true if chirality was applied successfully
      bool ApplyChirality
      (
        AtomVector< AtomComplete> &MOLECULE,
        const AtomVector< AtomComplete> &TEMPLATE,
        const bool &FORCE_ALTER_LABEL = false
      ) const;

      //! @brief Copy the double bond isometries of a different molecule onto this molecule
      //! This function assumes that MOLECULE currently has up-to-date double bond isometry information
      //! @param FORCE_ALTER_LABEL alter the chirality label even if the chirality cannot be adjusted via reflections
      //!        This is useful if later code has to assume that the mutation was successful
      //! @return true if isometry was applied successfully
      bool ApplyDoubleBondIsometry
      (
        AtomVector< AtomComplete> &MOLECULE,
        const AtomVector< AtomComplete> &TEMPLATE,
        const bool &FORCE_ALTER_LABEL = false
      ) const;

      //! @brief get the downstream atoms of a particular atom or bond
      //! @param ATOMA the first atom in the bond
      //! @param ATOMB the second atom in the bond
      //! @return everything that is reachable after going one way through bond A->B. Not valid for ring bonds
      const storage::Vector< size_t> &GetDownstreamAtoms( const size_t &ATOMA, const size_t &ATOMB) const;

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
    }; // class MutateChirality

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MUTATE_CHIRALITY_H_
