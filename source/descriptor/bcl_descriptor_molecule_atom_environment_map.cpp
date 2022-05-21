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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_molecule_atom_environment_map.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_environment.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "command/bcl_command_command_state.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // initialize static member data
    storage::HashMap< std::string, size_t> MoleculeAtomEnvironmentMap::s_AtomEnvironments = storage::HashMap< std::string, size_t>();
//    sched::Mutex& MoleculeAtomEnvironmentMap::GetMutex()
//    {
//      static sched::Mutex s_Mutex;
//      return s_Mutex;
//    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    MoleculeAtomEnvironmentMap::MoleculeAtomEnvironmentMap()
    :
       m_AtomEnvironmentsFilename( "chembl.h.atom_environment_hashmap.two_bonds.txt.gz"),
       m_RemoveH( false),
       m_EnvironmentSize( size_t( 2))
    {
    }

    //! @brief constructor
    MoleculeAtomEnvironmentMap::MoleculeAtomEnvironmentMap
    (
      const std::string &ENVIRONMENT_FILE,
      const bool &REMOVE_H,
      const size_t &ENVIRONMENT_SIZE
    )
    :
      m_AtomEnvironmentsFilename( ENVIRONMENT_FILE),
      m_RemoveH( REMOVE_H),
      m_EnvironmentSize( ENVIRONMENT_SIZE)
    {
      // Read in the atom environments file directly
      if( s_AtomEnvironments.IsEmpty())
      {
        io::IFStream( file);
        io::File::MustOpenIFStream
        (
          file,
          m_AtomEnvironmentsFilename
        );
        io::Serialize::Read( s_AtomEnvironments, file);
        io::File::CloseClearFStream( file);

        // Make sure the atom environemnts hashmap is not empty
        BCL_Assert( !s_AtomEnvironments.IsEmpty(), "Atom environment file is empty!");
      }
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeAtomEnvironmentMap
    MoleculeAtomEnvironmentMap *MoleculeAtomEnvironmentMap::Clone() const
    {
      return new MoleculeAtomEnvironmentMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeAtomEnvironmentMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeAtomEnvironmentMap::GetAlias() const
    {
      static const std::string s_name( "MoleculeAtomEnvironmentMap");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeAtomEnvironmentMap::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      // lock
//      GetMutex().Lock();

      util::SiPtr< const chemistry::ConformationInterface> mol_ptr( this->GetCurrentObject());
      chemistry::FragmentComplete mol( *mol_ptr);
      if( m_RemoveH)
      {
        mol.RemoveH();
      }

      // go over each atom in molecule
      for
      (
          auto itr_atoms( mol.GetAtomsIterator());
          itr_atoms.NotAtEnd();
          ++itr_atoms
      )
      {
        // generate atom environment string
        std::string atom_env_str;
        if( m_EnvironmentSize == size_t( 2))
        {
          atom_env_str = chemistry::AtomEnvironment::MakeAtomEnvironmentStringTwo( *itr_atoms);
        }
        else if( m_EnvironmentSize == size_t( 3))
        {
          atom_env_str = chemistry::AtomEnvironment::MakeAtomEnvironmentStringThree( *itr_atoms);
        }
        else if( m_EnvironmentSize == size_t( 4))
        {
          atom_env_str = chemistry::AtomEnvironment::MakeAtomEnvironmentStringFour( *itr_atoms);
        }
        else
        {
          BCL_Assert( true, "Invalid atom environemtn size specified. Please select either 2, 3, or 4");
        }

        // check if string is in reference hashmap
        if( s_AtomEnvironments.Find( atom_env_str) == s_AtomEnvironments.End())
        {
          // not found
          STORAGE( 0) = 0;
          break;
        }

        // all atom environments found
        STORAGE( 0) = 1;
      }
//      GetMutex().Unlock();
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeAtomEnvironmentMap::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // Load the default
      io::IFStream( file);
      io::File::MustOpenIFStream
      (
        file,
        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") +
        "atom_environments/" + ( m_AtomEnvironmentsFilename)
      );
      io::Serialize::Read( s_AtomEnvironments, file);
      io::File::CloseClearFStream( file);

      // Make sure the atom environments file actually contains some info
      BCL_Assert( !s_AtomEnvironments.IsEmpty(), "Atom environment file is empty!");
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeAtomEnvironmentMap::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Determines if molecule fragment radial atom environments exist in database"
      );
      parameters.AddInitializer
      (
        "exclude_hydrogen_atoms",
        "compare the local atom environments in the absence of hydrogen atoms; "
        "make sure to load a hashmap generated in the absence of explicit hydrogen atoms; "
        "faster but not guaranteed to perform as well",
        io::Serialization::GetAgent( &m_RemoveH),
        "false"
      );
      parameters.AddInitializer
      (
        "environment_size",
        "atoms within this number of bonds of an atom constitutes the local environment",
        io::Serialization::GetAgent( &m_EnvironmentSize),
        "2"
      );
      parameters.AddInitializer
      (
        "atom_environment_file",
        "the file containing the pre-generated fragment atom environment counts",
        io::Serialization::GetAgent( &m_AtomEnvironmentsFilename),
        m_RemoveH ?
        "chembl.noh.atom_environment_hashmap.two_bonds.txt.gz" :
        "chembl.h.atom_environment_hashmap.two_bonds.txt.gz"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
