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

#ifndef BCL_CHEMISTRY_MOLECULE_EVOLUTION_INFO_H_
#define BCL_CHEMISTRY_MOLECULE_EVOLUTION_INFO_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// headers from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
//#include "bcl_chemistry_molecule_evolutionary_optimizer.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_serializable_interface.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeEvolutionInfo
    //! @brief Convenience class that tracks pertinent information when evolving a molecule.
    //!
    //! @see @link example_chemistry_molecule_evolution_info.cpp @endlink
    //! @author geanesar, brownnp1
    //! @date Jun 18, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeEvolutionInfo :
      public util::SerializableInterface
      {

    private:

    /////////////
    // friends //
    /////////////

      friend class MoleculeEvolutionaryOptimizer; //!< Enable access to data to update molecules

    //////////
    // data //
    //////////

      // Name
      std::string m_Identifier;

      // Molecule whose evolutionary progress we are tracking
      FragmentComplete m_Molecule;

      // Fitness function value
      float m_Fitness;

      // Generational history
      storage::Vector< std::string> m_History;

      // Age in # generations
      size_t m_Age;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeEvolutionInfo();

      //! @brief constructor from data
      MoleculeEvolutionInfo
      (
        const std::string &IDENTIFIER,
        const FragmentComplete &MOLECULE,
        const float &FITNESS = util::GetUndefined< float>(),
        const storage::Vector< std::string> &HISTORY = storage::Vector< std::string>( 1, "Begin"),
        const size_t &AGE = 0
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeEvolutionInfo
      MoleculeEvolutionInfo *Clone() const
      {
        return new MoleculeEvolutionInfo( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of this class
      //! @return the name of this class
      const std::string &GetAlias() const;

      //! @brief return the molecule identifier
      //! @return the identifying string
      const std::string &GetMoleculeIdentifier() const;

      //! @brief get the stored molecule
      //! @return the member fragment;
      const FragmentComplete &GetMolecule() const;

    private:

      //! @brief get the stored molecule
      //! @return the member fragment;
      FragmentComplete GetMoleculeNonConst();

    public:

      //! @brief get the molecule fitness
      //! @return fitness value
      const float &GetMoleculeFitness() const;

      //! @brief get the molecule evolution history
      //! @return vector of strings dictating the history
      const storage::Vector< std::string> &GetMoleculeHistory() const;

      //! @brief get the generational age of the molecule
      //! @return molecule age
      const size_t &GetMoleculeAge() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief set the molecule identifier
      void SetMoleculeIdentifier( const std::string &IDENTIFIER);

      //! @brief set the stored molecule
      void SetMolecule( const FragmentComplete &MOL);

      //! @brief set the molecule fitness
      void SetMoleculeFitness( const float FITNESS);

      //! @brief set the molecule evolution history
      void SetMoleculeHistory( const storage::Vector< std::string> &HISTORY);

      //! @brief set the generational age of the molecule
      void SetMoleculeAge( const size_t AGE);

      //! @brief increment the molecule age by one generation
      void IncrementMoleculeAge();

      //! @brief append to molecule history by one generation
      void AppendToMoleculeHistory( const std::string &NEW_HISTORY_ENTRY);

    ///////////////
    // operators //
    ///////////////

      //! @brief less-than operator for MolInfos
      //! @return true if fitness of left operand is less than fitness of right operand
      bool operator <( const MoleculeEvolutionInfo &SECOND) const;

      //! @brief greater-than operator for MolInfos
      //! @return true if fitness of left operand is greater than fitness of right operand
      bool operator >( const MoleculeEvolutionInfo &SECOND) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief Set the members with LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      io::Serializer GetSerializer() const;

      }; // class MoleculeEvolutionInfo

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MOLECULE_EVOLUTION_INFO_H_
