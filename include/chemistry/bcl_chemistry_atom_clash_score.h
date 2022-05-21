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

#ifndef BCL_CHEMISTRY_ATOM_CLASH_SCORE_H_
#define BCL_CHEMISTRY_ATOM_CLASH_SCORE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_hash_map.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomClashScore
    //! @brief Scores clashes between atoms of a given conformation
    //! @details Potential that evalutes clashes in a given molecule conformation
    //!
    //! @see @link example_chemistry_atom_clash_score.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Aug 28, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomClashScore :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>,
      public signal::Slots
    {

    private:

      //! Whether to consider hydrogens
      bool m_ConsiderH;

      //! Whether to display the too-close atoms
      mutable bool m_Display;

      //! Molecule dependent variables. These variables are cached because they are expensive to compute and the same
      //! molecule may have the clash score called on it hundreds of times in succession
      mutable storage::Vector< AtomType>      m_LastAtomTypes;
      mutable storage::Vector< sdf::BondInfo> m_LastBondInfo;

      mutable float                           m_MaxVdwRadius;
      mutable linal::Matrix< float>           m_VdwSumMaxDistances;

      //! hash of molecules to preprocessed RMSD stuff
      mutable storage::HashMap
      <
        size_t,
        storage::Pair< double, storage::Vector< storage::Triplet< size_t, size_t, double> > >
      > m_Cache;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_HeavyInstance;

      //! @brief default constructor. Can take whether to consider hydrogens or tolerance
      //! @param HTOL hydrogen VDW-overlap tolerance. Default value causes non-covalent overlap of hydrogens to be ignored
      explicit AtomClashScore( const bool &CONSIDER_H = false, const bool &DISPLAY = false);

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AtomClashScore
      AtomClashScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetAlias() const;

      //! @brief set the tolerances and hydrogen consideration
      void Setup( const bool &CONSIDER_H)
      {
        if( m_ConsiderH != CONSIDER_H)
        {
          m_ConsiderH = CONSIDER_H;
          m_LastAtomTypes.Reset();
          m_LastBondInfo.Reset();
          m_Cache.Reset();
        }
      }

      void SetDisplay( const bool &DISP) const
      {
        m_Display = DISP;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate clashes for given atom pair
      //! @param MOLECULE molecule that needs to scored
      //! @return clash score for the given atom pair
      double operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief evaluate clashes for given atoms
      //! @param MOLECULE molecule that needs to scored
      //! @return clash score for the given atoms
      double operator()
      (
        const ConformationInterface &MOLECULE
      ) const;

      //! @brief get the minimum distance two particular atoms should be separated
      double GetMinimumNonClashingDistance( const size_t &A, const size_t &B) const
      {
        return m_VdwSumMaxDistances( A, B);
      }

      //! @brief get the clashing pairs
      const storage::Vector< storage::Triplet< size_t, size_t, double> > &GetClashingPairs
      (
        const ConformationInterface &MOLECULE
      ) const;

      //! @brief get the score for just a previous set of clashing pairs (saves refinding the clashing pairs, but may miss
      //!        some newly clashing pairs)
      storage::Vector< storage::Triplet< size_t, size_t, double> >::iterator UpdateClashingPairs
      (
        const AtomVector< AtomComplete> &MOLECULE,
        storage::Vector< storage::Triplet< size_t, size_t, double> > &CLASHING_PAIRS
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief Test whether this molecule is the same (constitutionally) as the molecule for which the state of this
      //!        class currently can handle
      bool IsMoleculeInfoCached( const ConformationInterface &CONF) const;

      //! @brief Update molecule change the molecule that this class will compute the clash score for
      //! @param MOL molecule of interest
      void UpdateMolecule( const ConformationInterface &CONF) const;

      //! @brief prepare the cache entry for a given molecule
      const storage::Pair< double, storage::Vector< storage::Triplet< size_t, size_t, double> > > &GetCachedInfo( const ConformationInterface &MOLECULE) const;

      //! @brief remove results for this object from the cache
      //! @param ARGUMENT argument that should be removed
      void RemoveResultFromCache( const ConformationInterface &ARGUMENT);

    }; // class AtomClashScore

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_CLASH_SCORE_H_
