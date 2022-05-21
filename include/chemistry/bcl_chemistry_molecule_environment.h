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

#ifndef BCL_CHEMISTRY_MOLECULE_ENVIRONMENT_H_
#define BCL_CHEMISTRY_MOLECULE_ENVIRONMENT_H_
// include the namespace headers
#include "bcl_chemistry.h"
#include "bcl_chemistry_conformation_comparison_interface.h"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_complete.h"
#include "bcl_chemistry_atom_environment_bender.h"
#include "bcl_chemistry_fragment_complete.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeEnvironment
    //! @brief stores the atoms whose distances from a atom of interest are no more than a certain number.
    //!
    //! @see @link example_chemistry_molecule_environment.cpp @endlink
    //! @author vuot2
    //! @date 06/29/2016
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MoleculeEnvironment :
      public ConformationComparisonInterface
    {

    public:
      typedef storage::Vector< AtomEnvironmentBender> t_MoleculeEnv;
      typedef AtomEnvironmentBender::AtomTypeEnum t_AtomTypeEnum;
      typedef AtomEnvironmentBender::Atom_type t_AtomType;
      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Elem;
      static const util::SiPtr< const util::ObjectInterface> s_ElemRC;
      static const util::SiPtr< const util::ObjectInterface> s_Atom;
      static const util::SiPtr< const util::ObjectInterface> s_AtomRC;

    private:
    //////////
    // data //
    //////////

      //! @brief atom type as string
      //! @param ATOM_TYPE the name of the atom type
      //! @return the string for the atom type
      static const std::string &GetAtomTypeName( const t_AtomType &ATOM_TYPE);
      t_MoleculeEnv m_MoleculeEnv;
      t_AtomTypeEnum m_AtomType;

    public:

    //////////////////
    // Constructions//
    //////////////////

      //! default constructor
      MoleculeEnvironment() {};

      //! constructor with Atom_type specification
      MoleculeEnvironment( const t_AtomType &ATOM_TYPE);

      //! @brief constructor for building an Molecule environment from bond distance limit, atom type, and fragment complete
      MoleculeEnvironment( const t_AtomType &ATOM_TYPE, const ConformationInterface &FRAGMENT);

      MoleculeEnvironment( const MoleculeEnvironment &OTHER);

      //! virtual copy constructor
      MoleculeEnvironment *Clone() const;

      //! @brief compares two atom environments
      bool operator ==( const MoleculeEnvironment &ATOM) const;

      //! @brief compares two atom environments
      bool operator !=( const MoleculeEnvironment &MOLECULE) const;

      //! @brief compute the tanimoto score between this and other molecule environment
      double operator()
      (
        const ConformationInterface& OTHER,
        const ConformationInterface& THIS
      ) const;

      //! instances of the class
      static const util::SiPtr< const util::ObjectInterface> s_Instances;

      //////////////////
      // data access ///
      //////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the atom type that this atom enviroment calculator calculates
      const t_AtomType &GetAtomType() const;

      //! @brief returns the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief creates atom environments from every atom of the molecule
      const t_MoleculeEnv &GetMoleculeEnvironment() const;

      // @brief Compute the TanimotoScore between this and OTHER
      double TanimotoScore( const MoleculeEnvironment &OTHER) const;

      // @brief Computes the Buser score between this and OTHER
      double BuserScore( const MoleculeEnvironment &OTHER) const;

    protected:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      // @brief converts vector of sorted objects into a map of objects and their counts
      storage::Map< AtomEnvironmentBender, size_t> ConvertVectorToMap() const;

      //! @brief adds the Key into MAP or increment its count( the value assiated with that key)
      static void AddToMap( const AtomEnvironmentBender &KEY, storage::Map< AtomEnvironmentBender, size_t> &MAP);
    };
  } // namespace chemistry
} // namespace bcl
#endif // BCL_CHEMISTRY_MOLECULE_ENVIRONMENT_H_
