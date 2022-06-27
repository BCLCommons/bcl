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
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "descriptor/bcl_descriptor_molecule_2da_smooth_sign_code.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_statistics.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_central_2da_sign.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "linal/bcl_linal_vector_operations.h"

using bcl::chemistry::AtomConformationalInterface;
using bcl::graph::Connectivity;

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Central2DASign::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Central2DASign()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Central2DASign::Central2DASign() :
      m_MaxBondDistanceFromCenter( size_t( 12)),
      m_Number2DASteps( size_t( 12)),
      m_Temperature( float( 100.0)),
      m_Smooth( true),
      m_2DASmoothSign( Molecule2DASmoothSignCode::e_2DASmoothSign)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {

        BCL_Assert
        (
          ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
          "Failed to create " + GetClassIdentifier()
        );
      }
    }

    //! @brief constructor from number of steps, and mapped atom property
    Central2DASign::Central2DASign
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t MAX_CENTER_BOND_DISTANCE,
      const size_t NUMBER_2DA_STEPS,
      const float &TEMPERATURE,
      const bool &SMOOTH
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_MaxBondDistanceFromCenter( MAX_CENTER_BOND_DISTANCE),
      m_Number2DASteps( NUMBER_2DA_STEPS),
      m_Temperature( TEMPERATURE),
      m_Smooth( SMOOTH),
      m_2DASmoothSign( Molecule2DASmoothSignCode::e_2DASmoothSign)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Central2DASign
    Central2DASign *Central2DASign::Clone() const
    {
      return new Central2DASign( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Central2DASign::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Central2DASign::GetAlias() const
    {
      static const std::string s_name( "Central2DASign");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Central2DASign::GetInternalDescriptors()
    {
      return
          iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
      (
        &m_2DASmoothSign,
        &m_2DASmoothSign + 1
      );
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Central2DASign::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      chemistry::FragmentComplete molecule( *current_mol);

      // Get one girth per atom
      storage::Vector< size_t> bond_girth( CalculateBondGirths( molecule));
      size_t min_bond_girth( math::Statistics::MinimumValue( bond_girth.Begin(), bond_girth.End()));

      //
      size_t i( 0);
      auto itr_desc( molecule.GetAtomsIterator());
      size_t inner_desc_size( m_2DASmoothSign.GetNormalSizeOfFeatures());
      for
      (
          auto itr_atoms( molecule.GetAtomsIterator());
          itr_atoms.NotAtEnd();
          ++itr_atoms, ++i, ++itr_desc
      )
      {
        size_t atoms_bond_girth( bond_girth( i));
        size_t atom_centrality( std::min( m_MaxBondDistanceFromCenter, size_t( atoms_bond_girth - min_bond_girth)));
        auto selected_spot( STORAGE.CreateSubVectorReference( inner_desc_size, inner_desc_size * atom_centrality));
        selected_spot += m_2DASmoothSign( itr_desc);
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return the bond girth from every atom
    //! @return bond girth at each atom index
    storage::Vector< size_t> Central2DASign::CalculateBondGirths( const chemistry::FragmentComplete &MOLECULE)
    {
      // get fragment graph
      if( m_Graph.GetSize() == size_t( 0))
      {
        chemistry::ConformationGraphConverter graph_maker;
        graph::ConstGraph< size_t, size_t> fragment_graph( graph_maker( MOLECULE));
        m_Graph = graph_maker( MOLECULE);
      }

      // get bond girth starting from each atom
      size_t mol_sz( MOLECULE.GetSize());
      storage::Vector< size_t> per_atom_bond_girths( mol_sz);
      for( size_t atom_i( 0); atom_i < mol_sz; ++atom_i)
      {
        // determine the distance from picked_atom vertex to all vertices
        util::ShPtr< storage::Vector< size_t> > bond_distances
        (
          graph::Connectivity::DistancesToOtherVertices( m_Graph, atom_i)
        );

        // keep the longest distance
        bond_distances->Sort( std::greater< size_t>());
        per_atom_bond_girths( atom_i) = bond_distances->FirstElement();
      }
      return per_atom_bond_girths;
    }

    //! @brief return the minimum bond girth
    //! @return shortest bond girth
    size_t Central2DASign::CalculateMinBondGirth( const chemistry::FragmentComplete &MOLECULE)
    {
      storage::Vector< size_t> bond_girths( CalculateBondGirths( MOLECULE));
      bond_girths.Sort( std::less< size_t>());
      return bond_girths( 0);
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Central2DASign::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // Setup 2DAs
      if( m_Smooth)
      {
        m_2DASmoothSign = Molecule2DASmoothSignCode( m_AtomProperty, m_Number2DASteps, m_Temperature, Molecule2DASmoothSignCode::e_2DASmoothSign);
      }
      else
      {
        m_2DASmoothSign = Molecule2DASmoothSignCode( m_AtomProperty, m_Number2DASteps, m_Temperature, Molecule2DASmoothSignCode::e_2DASign);
      }

      // Ensure bin size/resolution specified
      if( m_Number2DASteps == 0 || m_MaxBondDistanceFromCenter == 0.0)
      {
        ERR_STREAM << "Invalid step sizes specified";
        return false;
      }

      // Check that a property is provided
      if( m_AtomProperty.IsDefined() && m_AtomProperty->GetNormalSizeOfFeatures() != 1)
      {
        ERR_STREAM
        << "Expected a property that returned 1 properties per atom, but property returns "
        << m_AtomProperty->GetNormalSizeOfFeatures()
        << " values per atom ( property was "
        << m_AtomProperty->GetString()
        << ")";
        return false;
      }

      // Check that a property is provided
      if( !m_AtomProperty.IsDefined())
      {
        return false;
      }

      // return successfully if no issues
      return true;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Central2DASign::SetObjectHook()
    {
      m_Graph = graph::ConstGraph< size_t, size_t>();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Central2DASign::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Measure 2DASign of property at variable distances from molecule topological center");
      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the smooth radial distribution function",
        io::Serialization::GetAgent( &m_AtomProperty)
      );
      parameters.AddInitializer
      (
        "max_bonds_from_center",
        "size of outer array containing each signed 2DA",
        io::Serialization::GetAgentWithRange( &m_MaxBondDistanceFromCenter, 1, 1000000),
        "12"
      );
      parameters.AddInitializer
      (
        "2da_steps",
        "# of steps/bins (each of size = step size) used in the radial distribution function",
        io::Serialization::GetAgentWithRange( &m_Number2DASteps, 1, 1000000),
        "12"
      );
      parameters.AddInitializer
      (
        "smooth",
        "whether to smooth the 2DA",
        io::Serialization::GetAgent( &m_Smooth),
        "1"
      );
      if( m_Smooth)
      {
        parameters.AddInitializer
        (
          "temperature",
          "increasing temperature spreads autocorrelation across more distant bins",
          io::Serialization::GetAgentWithRange( &m_Temperature, 0.0, 1000.0),
          "100"
        );
      }
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
