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
#include "assemble/bcl_assemble_collector_aa_type.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorAAType::s_Instance
    (
      GetObjectInstances().AddInstance( new CollectorAAType())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorAAType::CollectorAAType() :
      m_AATypes( storage::Set< biol::AAType>::Create( biol::GetAATypes().ALA, biol::GetAATypes().ARG))
    {
    }

    //! @brief constructor taking member variables as paramters
    //! @param AA_TYPES set of aa types that will be collected from the list of residues
    CollectorAAType::CollectorAAType( const storage::Set< biol::AAType> &AA_TYPES) :
      m_AATypes( AA_TYPES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorAAType
    CollectorAAType *CollectorAAType::Clone() const
    {
      return new CollectorAAType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorAAType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &CollectorAAType::GetAlias() const
    {
      static const std::string s_Name( "CollectorAAType");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! Collect the residues that match the desired aa types from the provided residues
    //! @param RESIDUES entity that contains the residues to be collected from
    //! @return returns list of the collected residues that match the desired aa types
    util::SiPtrList< const biol::AABase>
    CollectorAAType::Collect( const util::SiPtrVector< const biol::AABase> &RESIDUES) const
    {
      // will hold all the collected residues that match the desired aa types
      util::SiPtrList< const biol::AABase> collected_aas;

      // iterate through the residues
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator aa_itr( RESIDUES.Begin()), aa_itr_end( RESIDUES.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // try to find the residue type in m_AATypes
        const storage::Set< biol::AAType>::const_iterator find_itr( m_AATypes.Find( ( *aa_itr)->GetType()));

        // true if the residue type was found in m_AATypes
        if( find_itr != m_AATypes.End())
        {
          // add the residue into collected_aas since it is one of the desired types
          collected_aas.PushBack( *aa_itr);
        }
      }

      // return the list of residues that were collected matching the desired types of residues
      return collected_aas;
    }

    //! @brief locate the residue that matches the desired aa type - if more than one just returns first encountered
    //! @param RESIDUES entity that contains the residues to be collected from
    //! @return const siptr to biol::AABase which is the first instance of the aa type found in RESIDUES
    util::SiPtr< const biol::AABase> CollectorAAType::Locate( const util::SiPtrVector< const biol::AABase> &RESIDUES) const
    {
      // get all the possible aas
      const util::SiPtrList< const biol::AABase> collected_aas( Collect( RESIDUES));

      // true if none were found
      if( collected_aas.IsEmpty())
      {
        // return empty siptr
        return util::SiPtr< const biol::AABase>();
      }

      // return first aa found
      return collected_aas.FirstElement();
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CollectorAAType::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AATypes, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CollectorAAType::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AATypes, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorAAType::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Collects residues which match any of the given aa types.");

      parameters.AddInitializer
      (
        "",
        "",
        io::Serialization::GetAgent( &m_AATypes),
        "(" + biol::GetAATypes().ALA->GetName() + ")"
      );

      return parameters;
    }

  } // namespace assemble

} // namespace bcl
