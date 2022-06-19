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
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_graph_marker.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_subgraph.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_ofstream.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_molecule_evolution_info.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_stopwatch.h"

// external includes

namespace bcl
{
  namespace chemistry
  {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    // brief constructor
    MoleculeEvolutionInfo::MoleculeEvolutionInfo() :
      m_Identifier(),
      m_Molecule(),
      m_Fitness(),
      m_History(),
      m_Age( 0)
    {
    }

    // brief constructor with arguments
    MoleculeEvolutionInfo::MoleculeEvolutionInfo
    (
      const std::string &IDENTIFIER,
      const FragmentComplete &MOLECULE,
      const float &FITNESS,
      const storage::Vector< std::string> &HISTORY,
      const size_t &AGE = 0
    ) :
      m_Identifier( IDENTIFIER),
      m_Molecule( MOLECULE),
      m_Fitness( FITNESS),
      m_History( HISTORY),
      m_Age( AGE)
    {
    }

    /////////////////
    // data access //
    /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoleculeEvolutionInfo::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of this class
    //! @return the name of this class
    const std::string &MoleculeEvolutionInfo::GetAlias() const
    {
      static const std::string s_name( "MoleculeEvolutionInfo");
      return s_name;
    }

    ////////////////
    // operations //
    ////////////////

    //! @brief less-than operator for MolInfos
    //! @return true if fitness of left operand is less than fitness of right operand
    bool MoleculeEvolutionInfo::operator <( const MoleculeEvolutionInfo &SECOND) const
    {
      return m_Fitness < SECOND.m_Fitness;
    }

    //! @brief greater-than operator for MolInfos
    //! @return true if fitness of left operand is greater than fitness of right operand
    bool MoleculeEvolutionInfo::operator >( const MoleculeEvolutionInfo &SECOND) const
    {
      return m_Fitness > SECOND.m_Fitness;
    }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // helper functions //
    //////////////////////

    //! @brief Set the members with LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeEvolutionInfo::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

    io::Serializer MoleculeEvolutionInfo::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Tracks the evolution of a single molecule.");

      return member_data;
    }

  } // namespace chemistry
} // namespace bcl
