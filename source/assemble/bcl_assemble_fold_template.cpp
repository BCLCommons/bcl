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
#include "assemble/bcl_assemble_fold_template.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "math/bcl_math_combination.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FoldTemplate::s_Instance
    (
      GetObjectInstances().AddInstance( new FoldTemplate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FoldTemplate::FoldTemplate() :
      m_Geometries(),
      m_PDBID( GetDefaultPDBID()),
      m_Topology(),
      m_TopologyCollector(),
      m_SortBySize( true)
    {
    }

    //! @brief construct from a vectors geometries
    //! @param GEOMETRIES bodies of SSEs in template
    //! @param TOPOLOGY_COLLECTOR collector to use to calculate the topology
    //! @param PDB_ID PDB ID of passed protein model
    //! @param SORT_BY_SIZE whether to sort the geometries by size
    FoldTemplate::FoldTemplate
    (
      const util::ShPtrVector< SSEGeometryPhiPsi> &GEOMETRIES,
      const util::ShPtr< CollectorTopologyInterface> &TOPOLOGY_COLLECTOR,
      const std::string &PDB_ID,
      const bool SORT_BY_SIZE
    ) :
      m_Geometries( GEOMETRIES),
      m_PDBID( PDB_ID),
      m_Topology(),
      m_TopologyCollector( TOPOLOGY_COLLECTOR),
      m_SortBySize( SORT_BY_SIZE)
    {
      // sort the geometries if requested
      if( m_SortBySize)
      {
        m_Geometries.Sort( SSEGeometryPhiPsiLessThan());
      }
    }

    //! @brief construct from a protein model and PDB ID, sorting the geometries by size
    //! @param PROTEIN_MODEL protein model
    //! @param TOPOLOGY_COLLECTOR collector to use to calculate the topology
    //! @param PDB_ID PDB ID of passed protein model
    //! @param TRIM_GEOMETRIES whether to trim the geometries by 1 residue at each terminus
    FoldTemplate::FoldTemplate
    (
      const ProteinModel &PROTEIN_MODEL,
      const util::ShPtr< CollectorTopologyInterface> &TOPOLOGY_COLLECTOR,
      const std::string &PDB_ID,
      const bool TRIM_GEOMETRIES
    ) :
      m_Geometries(),
      m_PDBID( PDB_ID),
      m_Topology(),
      m_TopologyCollector( TOPOLOGY_COLLECTOR),
      m_SortBySize( true)
    {
      // create the geometries from the protein SSEs
      CreateGeometries( PROTEIN_MODEL, TRIM_GEOMETRIES);
    }

    //! @brief copy constructor
    //! @param FOLD_TEMPLATE FoldTemplate to be copied
    FoldTemplate::FoldTemplate( const FoldTemplate &FOLD_TEMPLATE) :
      m_Geometries( FOLD_TEMPLATE.m_Geometries.HardCopy()),
      m_PDBID( FOLD_TEMPLATE.m_PDBID),
      // m_topology is not copied because it will still have pointers to FOLD_TEMPLATE
      // it will be recalculated the first time it is needed
      m_Topology(),
      m_TopologyCollector( FOLD_TEMPLATE.m_TopologyCollector),
      m_SortBySize( FOLD_TEMPLATE.m_SortBySize)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldTemplate
    FoldTemplate *FoldTemplate::Clone() const
    {
      return new FoldTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FoldTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief calls GetIndentification on each geometry
    //! @return std::string of GetIndentification from each geometry, separated by a new line
    std::string FoldTemplate::GetGeometryIndentifications() const
    {
      // initialize return string
      std::string geometry_string;

      // iterate through the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geo_itr( m_Geometries.Begin()),
          geo_itr_end( m_Geometries.End());
        geo_itr != geo_itr_end;
        ++geo_itr
      )
      {
        // add the identification to the string
        geometry_string += ( *geo_itr)->GetIdentification() + "\n";
      }

      // end
      return geometry_string;
    }

    //! @brief get the helical geometries
    //! @return the helical geometries
    util::SiPtrVector< const SSEGeometryPhiPsi> FoldTemplate::GetHelicalGeometries() const
    {
      // initialize geometry vector
      util::SiPtrVector< const SSEGeometryPhiPsi> helical_geometries;

      // iterate through the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geo_itr( m_Geometries.Begin()),
          geo_itr_end( m_Geometries.End());
        geo_itr != geo_itr_end;
        ++geo_itr
      )
      {
        // if the geometry is a helix
        if( ( *geo_itr)->GetType() == biol::GetSSTypes().HELIX)
        {
          // add to the geometry vector
          helical_geometries.PushBack( **geo_itr);
        }
      }

      // end
      return helical_geometries;
    }

    //! @brief get the strand geometries
    //! @return the strand geometries
    util::SiPtrVector< const SSEGeometryPhiPsi> FoldTemplate::GetStrandGeometries() const
    {
      // initialize geometry vector
      util::SiPtrVector< const SSEGeometryPhiPsi> strand_geometries;

      // iterate through the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geo_itr( m_Geometries.Begin()),
          geo_itr_end( m_Geometries.End());
        geo_itr != geo_itr_end;
        ++geo_itr
      )
      {
        // if the geometry is a strand
        if( ( *geo_itr)->GetType() == biol::GetSSTypes().STRAND)
        {
          // add to the geometry vector
          strand_geometries.PushBack( **geo_itr);
        }
      }

      // end
      return strand_geometries;
    }

    //! @brief returns whether the topology is initialized or not
    //! @return whether the topology is initialized or not
    bool FoldTemplate::IsTopologyInitialized() const
    {
      return ( m_Topology.GetGraph().GetNumberVertices() > 0);
    }

    //! @brief creates the topology from the geometry information
    //! @return the topology
    void FoldTemplate::CalculateTopology()
    {
      // calculate the SSE geometries using the phi/psi angles if needed
      if( !HasDefinedGeometries())
      {
        CalculateSSEGeometries();
      }

      // then initialize the topology graph
      m_Topology = m_TopologyCollector->CalculateTopology( m_Geometries);
    }

    //! @brief set the topology
    //! @param TOPOLOGY Topology to be copied
    void FoldTemplate::SetTopology( const Topology &TOPOLOGY)
    {
      m_Topology = TOPOLOGY;
    }

    //! @brief get the center of mass of the fold template
    //! @return the center of mass of the fold template
    linal::Vector3D FoldTemplate::GetCenter() const
    {
      // get the coordinates
      storage::Vector< linal::Vector3D> coords( CollectCoordinates());

      // calculate center of mass and return
      return coord::CenterOfMass( util::ConvertToSiPtrVector( coords));
    }

    //! @brief get the radius of gyration of the fold template
    //! @return the radius of gyration of the fold template
    double FoldTemplate::GetRadiusOfGyration() const
    {
      // get the coordinates
      storage::Vector< linal::Vector3D> coords( CollectCoordinates());

      // calculate the radius of gyration and return
      return coord::RadiusOfGyration( util::ConvertToSiPtrVector( coords));
    }

    //! @brief gets a subdomain of the appropriate size
    //! @param HELICES number of helices in the subdomain
    //! @param STRANDS number of strands in the subdomain
    //! @return a subdomain of the appropriate size
    FoldTemplate FoldTemplate::GetSubDomain( const size_t HELICES, const size_t STRANDS) const
    {
      // get all the possible subdomains
      util::ShPtrList< FoldTemplate> all_subdomains( GetAllSubdomains( HELICES, STRANDS));

      // exit with an empty fold template if no subdomains could be identified
      if( all_subdomains.IsEmpty())
      {
        BCL_MessageStd( "Unable to find a subdomain with the requested SSE composition");
        return FoldTemplate();
      }

      // get a random subdomain
      FoldTemplate random_template
      (
        **random::GetGlobalRandom().Iterator( all_subdomains.Begin(), all_subdomains.End(), all_subdomains.GetSize())
      );

      // calculate the topology
      random_template.CalculateTopology();

      // end
      return random_template;
    }

    //! @brief gets a subdomain that matches the passed SSEs by given comparison
    //! @param SSES sses used to chose subdomains based on length
    //! @param SSE_GEOMETRY_COMPARE comparison method
    //! @return a subdomain that matches the passed SSEs by given comparison
    FoldTemplate FoldTemplate::GetSubDomain
    (
      const util::SiPtrVector< const SSE> &SSES,
      const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE
    ) const
    {
      // get the helical sses
      const CollectorSSE helix_collector( biol::GetSSTypes().HELIX);
      const util::SiPtrList< const SSE> helices( helix_collector.Collect( SSES));

      // get the strand sses
      const CollectorSSE strand_collector( biol::GetSSTypes().STRAND);
      const util::SiPtrList< const SSE> strands( strand_collector.Collect( SSES));

      // get all the subdomains of the right SSE composition
      util::ShPtrList< FoldTemplate> all_subdomains( GetAllSubdomains( helices.GetSize(), strands.GetSize()));

      // exit with an empty fold template if no subdomains could be identified
      if( all_subdomains.IsEmpty())
      {
        BCL_MessageVrb
        (
          "Unable to find a subdomain with the requested SSE composition: H(" + util::Format()( helices.GetSize()) +
          ") S(" + util::Format()( strands.GetSize()) + ")"
        );
        return FoldTemplate();
      }

      // initialize probability distribution to hold templates that satisfy the SSE lengths, they are weighted by the
      // square of the average neighbor count
      util::ShPtrVector< FoldTemplate> suitable_subdomains;

      // iterate through the list of all subdomains
      for
      (
        util::ShPtrList< FoldTemplate>::iterator template_itr( all_subdomains.Begin()),
          template_itr_end( all_subdomains.End());
        template_itr != template_itr_end; ++template_itr
      )
      {
        // if the geometry sizes match the sse sizes
        if( ( *template_itr)->HasSimilarSizeGeometries( SSES, SSE_GEOMETRY_COMPARE))
        {
          // add this subdomain
          suitable_subdomains.PushBack( *template_itr);
        }
      }

      // exit with an empty fold template if no subdomains could be identified
      if( suitable_subdomains.IsEmpty())
      {
        BCL_MessageVrb
        (
          "Unable to find a subdomain with the requested SSE composition and sizes"
        );
        return FoldTemplate();
      }

      // get a random subdomain
      FoldTemplate random_template
      (
        **random::GetGlobalRandom().Iterator
        (
          suitable_subdomains.Begin(), suitable_subdomains.End(), suitable_subdomains.GetSize()
        )
      );

      // calculate the topology
      random_template.CalculateTopology();

      // end
      return random_template;
    }

    //! @brief gets the default PDB ID to be used if none is supplied
    //! @return the default PDB ID to be used if none is supplied
    const std::string &FoldTemplate::GetDefaultPDBID()
    {
      static const std::string s_default_pdb_id( "DFLT");

      return s_default_pdb_id;
    }

    //! @brief returns whether the geometries are close in size to the passed SSEs
    //! @param SSES sses to compare geometries to
    //! @param GEOMETRY_COMPARE function to decide whether an SSE should be fit into a geometry
    //! @return whether the geometries are close in size to the passed SSEs
    bool FoldTemplate::HasSimilarSizeGeometries
    (
      const util::SiPtrVector< const SSE> &SSES,
      const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &GEOMETRY_COMPARE
    ) const
    {
      // get the helical sses
      const CollectorSSE helix_collector( biol::GetSSTypes().HELIX);
      util::SiPtrList< const SSE> helices( helix_collector.Collect( SSES));

      // get the strand sses
      const CollectorSSE strand_collector( biol::GetSSTypes().STRAND);
      util::SiPtrList< const SSE> strands( strand_collector.Collect( SSES));

      // get the geometries from this fold template
      const util::SiPtrVector< const SSEGeometryPhiPsi> helical_geometries( GetHelicalGeometries());
      const util::SiPtrVector< const SSEGeometryPhiPsi> strand_geometries( GetStrandGeometries());

      // if the number of SSEs and geometries don't match
      if( helices.GetSize() != helical_geometries.GetSize() || strands.GetSize() != strand_geometries.GetSize())
      {
        return false;
      }

      helices.Sort( SSELessThanBySize());
      strands.Sort( SSELessThanBySize());

      // get iterators on the helical geometries
      util::SiPtrVector< const SSEGeometryPhiPsi>::const_iterator geometry_itr( helical_geometries.Begin());
      const util::SiPtrVector< const SSEGeometryPhiPsi>::const_iterator geometry_itr_end( helical_geometries.End());

      // initialize bool that wrong size was found
      bool found_wrong_size( false);

      // iterate through the helical SSEs and geometries
      for
      (
        util::SiPtrList< const SSE>::const_iterator sse_itr( helices.Begin()), sse_itr_end( helices.End());
        sse_itr != sse_itr_end && geometry_itr != geometry_itr_end; ++sse_itr, ++geometry_itr
      )
      {
        // if the compare function returns false
        if( !GEOMETRY_COMPARE( **sse_itr, **geometry_itr))
        {
          // set the bool and break
          found_wrong_size = true;
          break;
        }
      }

      // only look at the strand geometries if the helix geometries matched
      if( !found_wrong_size)
      {
        // get iterators on the strand geometries
        util::SiPtrVector< const SSEGeometryPhiPsi>::const_iterator geometry_itr( strand_geometries.Begin());
        const util::SiPtrVector< const SSEGeometryPhiPsi>::const_iterator geometry_itr_end( strand_geometries.End());

        // iterate through the strand SSEs and geometries
        for
        (
          util::SiPtrList< const SSE>::const_iterator sse_itr( strands.Begin()), sse_itr_end( strands.End());
          sse_itr != sse_itr_end && geometry_itr != geometry_itr_end; ++sse_itr, ++geometry_itr
        )
        {
          // if the compare function returns false
          if
          (
            !GEOMETRY_COMPARE( **sse_itr, **geometry_itr)
          )
          {
            // set the bool and break
            found_wrong_size = true;
            break;
          }
        }
      }

      // end
      return !found_wrong_size;
    }

    //! @brief returns whether there is at least one geometry that is similar to the SSE
    //! @param SP_SSE sse to compare geometries to
    //! @param GEOMETRY_COMPARE function to decide whether an SSE should be fit into a geometry
    //! @return whether there is at least one geometry that is similar to the SSE
    bool FoldTemplate::HasSimilarSizeGeometry
    (
      const util::SiPtr< const SSE> &SP_SSE,
      const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &GEOMETRY_COMPARE
    ) const
    {
      // get the SSE type
      const biol::SSType ss_type( SP_SSE->GetType());

      // if it is not helix or strand return fals
      if( !ss_type.IsDefined())
      {
        return false;
      }

      // get the appropriate geometries
      const util::SiPtrVector< const SSEGeometryPhiPsi> template_geometries
      (
        ss_type == biol::GetSSTypes().HELIX ? GetHelicalGeometries() : GetStrandGeometries()
      );

      // iterate over the geometries
      for
      (
          util::SiPtrVector< const SSEGeometryPhiPsi>::const_iterator geometry_itr( template_geometries.Begin()),
          geometry_itr_end( template_geometries.End());
        geometry_itr != geometry_itr_end; ++geometry_itr
      )
      {
        // if the geometry matches
        if( GEOMETRY_COMPARE( *SP_SSE, **geometry_itr))
        {
          return true;
        }
      }

      // if this point is reached no matching geometries were found so return false
      return false;
    }

    //! @brief returns whether all the geometries are defined
    //! @return whether all the geometries are defined
    bool FoldTemplate::HasDefinedGeometries() const
    {
      // iterate over the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geometry_itr( m_Geometries.Begin()),
          geometry_itr_end( m_Geometries.End());
        geometry_itr != geometry_itr_end; ++geometry_itr
      )
      {
        // if the geometry is undefined
        if( !( *geometry_itr)->IsDefined())
        {
          return false;
        }
      }

      // all geometries are defined
      return true;
    }

    //! @brief returns whether all the angles are defined
    //! @return whether all the angles are defined
    bool FoldTemplate::HasDefinedAngles() const
    {
      // iterate over the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geometry_itr( m_Geometries.Begin()),
          geometry_itr_end( m_Geometries.End());
        geometry_itr != geometry_itr_end; ++geometry_itr
      )
      {
        // if the angle is undefined
        if( !( *geometry_itr)->GetPhiPsi()->IsDefined())
        {
          return false;
        }
      }

      // all angles are defined
      return true;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief performs translation given by linal::Vector3D
    //! @param TRANSLATION linal::Vector3D to be used for translation
    void FoldTemplate::Translate( const linal::Vector3D &TRANSLATION)
    {
      // create a transformation matrix from the translation
      math::TransformationMatrix3D transform;
      transform( TRANSLATION);

      // apply translation
      Transform( transform);
    }

    //! @brief performs rotation & translation given by 4D transformation matrix
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be used for transformation
    void FoldTemplate::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      // iterate through the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::iterator geo_itr( m_Geometries.Begin()),
          geo_itr_end( m_Geometries.End());
        geo_itr != geo_itr_end;
        ++geo_itr
      )
      {
        // apply transformation
        ( *geo_itr)->Transform( TRANSFORMATION_MATRIX_3D);
      }
    }

    //! @brief performs rotation given by math::RotationMatrix3D
    //! @param ROTATIONMATRIX3D math::RotationMatrix3D to be used for rotation
    void FoldTemplate::Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D)
    {
      // apply rotation
      Transform( math::TransformationMatrix3D( ROTATIONMATRIX3D));
    }

    //! @brief fits the given SSEs into this fold template and returns a domain
    //! @param SSES_TO_FIT vector of SSEs to place into the fold template
    //! @return a domain with the placed SSES
    Domain FoldTemplate::FitSSEs( const util::SiPtrVector< const SSE> &SSES_TO_FIT) const
    {
      // initialize vectors of helices and strands
      util::ShPtrVector< SSE> helices;
      util::ShPtrVector< SSE> strands;

      // iterate through the SSEs
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( SSES_TO_FIT.Begin()), sse_itr_end( SSES_TO_FIT.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // add to the appropriate vector
        if( ( *sse_itr)->GetType() == biol::GetSSTypes().HELIX)
        {
          helices.PushBack( util::ShPtr< SSE>( ( *sse_itr)->Clone()));
        }
        else if( ( *sse_itr)->GetType() == biol::GetSSTypes().STRAND)
        {
          strands.PushBack( util::ShPtr< SSE>( ( *sse_itr)->Clone()));
        }
      }

      // sort by size
      std::sort( helices.Begin(), helices.End(), SSELessThanBySize());
      std::sort( strands.Begin(), strands.End(), SSELessThanBySize());

      // get the helix and strand geometries
      util::SiPtrVector< const SSEGeometryPhiPsi> helical_geometries( GetHelicalGeometries());
      util::SiPtrVector< const SSEGeometryPhiPsi> strand_geometries( GetStrandGeometries());

      // make sure the vectors are the same size
      MatchSSEsAndGeometries( helices, helical_geometries);
      MatchSSEsAndGeometries( strands, strand_geometries);

      // create iterators on the SSEs
      util::ShPtrVector< SSE>::iterator helix_itr( helices.Begin());
      const util::ShPtrVector< SSE>::iterator helix_itr_end( helices.End());
      util::ShPtrVector< SSE>::iterator strand_itr( strands.Begin());
      const util::ShPtrVector< SSE>::iterator strand_itr_end( strands.End());

      // create new domain sse vector
      util::ShPtrVector< SSE> domain_sses;

      // iterate through the helices
      for
      (
        util::SiPtrVector< const SSEGeometryPhiPsi>::const_iterator geometry_itr( helical_geometries.Begin()),
          geometry_itr_end( helical_geometries.End());
        geometry_itr != geometry_itr_end && helix_itr != helix_itr_end;
        ++geometry_itr, ++helix_itr
      )
      {
        // apply the fit
        biol::AASequenceFactory::FitSequence( **helix_itr, *( ( *geometry_itr)->GetPhiPsi()), biol::GetSSTypes().HELIX);

        // update the SSE geometry since the associated sequence has changed
        ( *helix_itr)->SetGeometry();

        // add the sse to the domain vector
        domain_sses.PushBack( *helix_itr);
      }

      // iterate through the strands
      for
      (
        util::SiPtrVector< const SSEGeometryPhiPsi>::const_iterator geometry_itr( strand_geometries.Begin()),
          geometry_itr_end( strand_geometries.End());
        geometry_itr != geometry_itr_end && strand_itr != strand_itr_end;
        ++geometry_itr, ++strand_itr
      )
      {
        // apply the fit
        biol::AASequenceFactory::FitSequence
        (
          **strand_itr,
          *( ( *geometry_itr)->GetPhiPsi()),
          biol::GetSSTypes().STRAND
        );

        // update the SSE geometry since the associated sequence has changed
        ( *strand_itr)->SetGeometry();

        // add the sse to the domain vector
        domain_sses.PushBack( *strand_itr);
      }

      // end
      return Domain( domain_sses);
    }

    //! @brief calculates the SSEGeometries using the associated phi/psi information
    void FoldTemplate::CalculateSSEGeometries()
    {
      // iterate through the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::iterator geo_itr( m_Geometries.Begin()),
          geo_itr_end( m_Geometries.End());
        geo_itr != geo_itr_end; ++geo_itr
      )
      {
        // set it using the phi/psi angles
        ( *geo_itr)->SetSSEGeometryUsingPhiPsi();
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param FOLD_TEMPLATE FoldTemplate to be copied
    //! @return This template after all members are assigned to values from FOLD_TEMPLATE
    FoldTemplate &FoldTemplate::operator =( const FoldTemplate &FOLD_TEMPLATE)
    {
      // m_topology is not copied because it will still have pointers to FOLD_TEMPLATE
      // it will be recalculated the first time it is needed
      m_Geometries = FOLD_TEMPLATE.m_Geometries.HardCopy();
      m_PDBID = FOLD_TEMPLATE.m_PDBID;
      m_TopologyCollector = FOLD_TEMPLATE.m_TopologyCollector;
      m_SortBySize = FOLD_TEMPLATE.m_SortBySize;

      // return updated template
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FoldTemplate::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Geometries, ISTREAM);
      io::Serialize::Read( m_PDBID, ISTREAM);
      io::Serialize::Read( m_TopologyCollector, ISTREAM);
      io::Serialize::Read( m_SortBySize, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FoldTemplate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Geometries, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PDBID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TopologyCollector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SortBySize, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write minimal information needed to construct a FoldTemplate to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @return outputstream which was written to
    std::ostream &FoldTemplate::WriteCompact( std::ostream &OSTREAM) const
    {
      // write out the PDB ID number of geometries
      OSTREAM << m_PDBID << ' ' << util::Format()( m_Geometries.GetSize()) << '\n';

      // write out just the phi/psi information by iterating over the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geo_itr( m_Geometries.Begin()),
          geo_itr_end( m_Geometries.End());
        geo_itr != geo_itr_end; ++geo_itr
      )
      {
        // write out the type, coords, and size
        OSTREAM << ( *geo_itr)->GetType().GetName() << ' '
                << WriteCoordinates( ( *geo_itr)->GetPhiPsi()->GetN())
                << WriteCoordinates( ( *geo_itr)->GetPhiPsi()->GetCA())
                << WriteCoordinates( ( *geo_itr)->GetPhiPsi()->GetC())
                << util::Format()( ( *geo_itr)->GetPhiPsi()->GetAngles().GetSize());

        // get iterators on the phi psi angles
        storage::Vector< storage::VectorND< 2, double> >::const_iterator angle_itr
        (
          ( *geo_itr)->GetPhiPsi()->GetAngles().Begin()
        );
        storage::Vector< storage::VectorND< 2, double> >::const_iterator angle_itr_end
        (
          ( *geo_itr)->GetPhiPsi()->GetAngles().End()
        );

        // write out the first psi (phi is nan)
        OSTREAM << ' ' << util::Format()( angle_itr->Second());
        ++angle_itr;
        --angle_itr_end;

        // iterate over the phi/psi angles
        for( ; angle_itr != angle_itr_end; ++angle_itr)
        {
          // write the phi and psi angles
          OSTREAM << ' ' << util::Format()( angle_itr->First()) << ' ' << util::Format()( angle_itr->Second());
        }

        // write out the last phi (psi is nan)
        OSTREAM << ' ' << util::Format()( ( *geo_itr)->GetPhiPsi()->GetAngles().LastElement().First());

        // move to next line
        OSTREAM << '\n';
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create geometry member variables
    //! @param PROTEIN_MODEL protein model
    //! @param TRIM_GEOMETRIES whether to trim the geometries
    void FoldTemplate::CreateGeometries( const ProteinModel &PROTEIN_MODEL, const bool TRIM_GEOMETRIES)
    {
      // store SSEs
      util::SiPtrVector< const SSE> sses_vector
      (
        PROTEIN_MODEL.GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
      );

      // if sorting is selected
      if( m_SortBySize)
      {
        std::sort( sses_vector.Begin(), sses_vector.End(), SSELessThanBySize());
      }

      // initialize a geometries vector
      util::ShPtrVector< SSEGeometryPhiPsi> geometries;

      // iterate through SSEs in PROTEIN_MODEL
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sses_vector.Begin()), sse_itr_end( sses_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if the geometries should be trimmed
        if( TRIM_GEOMETRIES)
        {
          // make a copy of the SSE
          SSE copy_sse( **sse_itr);

          // clip the ends
          copy_sse.ClipEnds( 1);

          // update the geometry since ClipEnds is an AASequence function
          copy_sse.SetGeometry();

          // add to the geometries
          geometries.PushBack( util::ShPtr< SSEGeometryPhiPsi>( new SSEGeometryPhiPsi( copy_sse)));
        }
        else
        {
          // add the original SSE to the geometries
          geometries.PushBack( util::ShPtr< SSEGeometryPhiPsi>( new SSEGeometryPhiPsi( **sse_itr)));
        }
      }

      // initialize the topology
      // then initialize the topology graph
      m_Topology = m_TopologyCollector->CalculateTopology( geometries);

      // iterate over the current geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geo_itr( geometries.Begin()),
          geo_itr_end( geometries.End());
        geo_itr != geo_itr_end; ++geo_itr
      )
      {
        // if this geometry does also exist in the topology as a vertex
        // meaning it has at least one SSEPacking that is valid according to criteria and therefore has an edge
        if( m_Topology.GetGraph().HasVertex( *geo_itr))
        {
          m_Geometries.PushBack( *geo_itr);
        }
        else
        {
          BCL_MessageStd( util::Format()( ( *geo_itr)->GetIdentification()) + " has no neighbors");
        }
      }
    }

    //! @brief calculate coordinates for each body in a fold template
    //! @return coordinates for a fold template
    storage::Vector< linal::Vector3D> FoldTemplate::CollectCoordinates() const
    {
      // initialize vector of coordinates
      storage::Vector< linal::Vector3D> coords;

      // iterate through the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geo_itr( m_Geometries.Begin()),
          geo_itr_end( m_Geometries.End());
        geo_itr != geo_itr_end;
        ++geo_itr
      )
      {
        // get the fragment geometries
        util::SiPtrVector< const coord::GeometryInterface> fragments( ( *geo_itr)->GetGeometries());

        // if the fragment list is empty
        if( fragments.IsEmpty())
        {
          // pushback the main geometry
          fragments.PushBack( **geo_itr);
        }

        // iterate through the fragments
        for
        (
          util::SiPtrVector< const coord::GeometryInterface>::const_iterator fragment_itr( fragments.Begin()),
            fragment_itr_end( fragments.End());
          fragment_itr != fragment_itr_end;
          ++fragment_itr
        )
        {
          // add the begin, center, and end points to the coordinates
          coords.PushBack( ( *fragment_itr)->GetMainAxis().GetStartPoint());
          coords.PushBack( ( *fragment_itr)->GetCenter());
          coords.PushBack( ( *fragment_itr)->GetMainAxis().GetEndPoint());
        }
      }

      // end
      return coords;
    }

    //! @brief writes the given coordinates to one line
    //! @param COORDS Vector3D to be written
    //! @return string containing the coordinates on one line
    std::string FoldTemplate::WriteCoordinates( const linal::Vector3D &COORDS)
    {
      // initialize the string
      std::string coords_string;

      // write the orientation
      for( size_t i( 0); i != 3; ++i)
      {
        coords_string += util::Format()( COORDS( i)) + " ";
      }

      // end
      return coords_string;
    }

    //! @brief trims vectors of SSEs and geometries (whichever is biggest so that they match)
    //! @param SSES vector of SSEs
    //! @param GEOS vector of SSEGeometries
    void FoldTemplate::MatchSSEsAndGeometries
    (
      util::ShPtrVector< SSE> &SSES, util::SiPtrVector< const SSEGeometryPhiPsi> &GEOS
    )
    {
      // get the sizes of the vectors
      const size_t nr_sses( SSES.GetSize());
      const size_t nr_geos( GEOS.GetSize());

      // there are more SSEs than geometries
      if( nr_sses > nr_geos)
      {
        // iterate over the size difference
        for( size_t i( 0), i_end( nr_sses - nr_geos); i != i_end; ++i)
        {
          // remove a random element
          SSES.RemoveRandomElement();
        }
      }
      // there are more geometries than SSEs
      else if( nr_sses < nr_geos)
      {
        // iterate over the size difference
        for( size_t i( 0), i_end( nr_geos - nr_sses); i != i_end; ++i)
        {
          // remove a random element
          GEOS.RemoveRandomElement();
        }
      }
    }

    //! @brief get the average number of close geometries
    //! @return the average number of close geometries
    double FoldTemplate::AverageNumberOfCloseGeometries() const
    {
      // initialize sum
      double neighbor_sum( 0.0);

      // iterate through the geometries
      for
      (
        util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geometry_itr( m_Geometries.Begin()),
          geometry_itr_end( m_Geometries.End());
        geometry_itr != geometry_itr_end; ++geometry_itr
      )
      {
        // add to the number of neighbors
        neighbor_sum += NumberOfCloseGeometries( *geometry_itr);
      }

      // get the number of geometries
      const double nr_geometries( m_Geometries.GetSize());

      // end
      return nr_geometries == 0.0 ? 0.0 : neighbor_sum / nr_geometries;
    }

    //! @brief get the number of geometries within DISTANCE
    //! @param GEOMETRY geometry of interest
    //! @return the number of bodies within DISTANCE
    size_t FoldTemplate::NumberOfCloseGeometries( const util::ShPtr< SSEGeometryPhiPsi> &GEOMETRY) const
    {
      return GetNeighborGeometries( GEOMETRY).GetSize();
    }

    //! @brief get geometries that contact the passed geometry in the fold template
    //! @param GEOMETRY geometry of interest
    //! @return ShPtrVector of geometries that contact the passed geometry in the fold template
    util::ShPtrVector< SSEGeometryPhiPsi> FoldTemplate::GetNeighborGeometries
    (
      const util::ShPtr< SSEGeometryPhiPsi> &GEOMETRY
    ) const
    {
      // initialize geometry vector
      util::ShPtrVector< SSEGeometryPhiPsi> neighbors;

      // get the corresponding vertex from the topology graph
      util::ShPtr< Topology::GraphType::VertexType> vertex
      (
        m_Topology.GetGraph().FindVertex( util::SiPtr< const SSEGeometryInterface>( GEOMETRY))
      );

      // is the ShPtr is undefined
      if( !vertex.IsDefined())
      {
        // return an empty vector
        BCL_MessageStd
        (
          "SSEGeometry " + GEOMETRY->GetIdentification() +
          " was not found in the fold template topology, returning empty neighbor geometry vector"
        );

        return neighbors;
      }

      // iterate through the edges
      for
      (
        Topology::GraphType::EdgeContainerType::const_iterator edge_itr( vertex->GetEdges().Begin()),
          edge_itr_end( vertex->GetEdges().End());
        edge_itr != edge_itr_end; ++edge_itr
      )
      {
        // get a pointer to the target data
        util::SiPtr< const SSEGeometryInterface> target_geometry( edge_itr->GetTarget()->GetData());

        // iterate through the geometries
        for
        (
          util::ShPtrVector< SSEGeometryPhiPsi>::const_iterator geometry_itr( m_Geometries.Begin()),
            geometry_itr_end( m_Geometries.End());
          geometry_itr != geometry_itr_end; ++geometry_itr
        )
        {
          // check if the geometry pointer is also pointing to the target data
          if( geometry_itr->GetPointer() == target_geometry.GetPointer())
          {
            // add to the neighbor vector and continue
            neighbors.PushBack( *geometry_itr);
            continue;
          }
        }
      }

      // end
      return neighbors;
    }

    //! @brief gets all the possible subdomains of the requested size
    //! @param HELICES number of helices in the subdomain
    //! @param STRANDS number of strands in the subdomain
    //! @return all the possible subdomains of the requested size
    util::ShPtrList< FoldTemplate> FoldTemplate::GetAllSubdomains( const size_t HELICES, const size_t STRANDS) const
    {
      // check if the requested subdomain size is larger than the fold template
      if( HELICES > GetHelicalGeometries().GetSize() || STRANDS > GetStrandGeometries().GetSize())
      {
        // warn the user and return an empty template
        BCL_MessageVrb
        (
          "Unable to find a subdomain of size " + util::Format()( HELICES) +
          " and " + util::Format()( STRANDS) + " because the fold template is too small"
        );
        return util::ShPtrList< FoldTemplate>();
      }

      // The next step is to find a subgraph from the graph that has the requested number of helices and strands.
      // The approach here is to first store all possible combinations of helices and strands that can be generated to
      // match the number requested (HELICES & STRANDS).  Then each strand combination is added to each helix
      // combination to get HELICES * STRANDS total collections of SSEs.  Next, a subgraph is generated for each
      // collection.  Only connected subgraphs (meaning each SSE is in contact with at least one other SSE) is retained.
      // A fold template is then generated from one of these subgraphs randomly.

      // get the helix and strand geometries
      util::SiPtrVector< const SSEGeometryPhiPsi> helix_geometries( GetHelicalGeometries());
      util::SiPtrVector< const SSEGeometryPhiPsi> strand_geometries( GetStrandGeometries());

      // construct combinations for the helix and strand geometries
      math::Combination< util::SiPtr< const SSEGeometryPhiPsi> > helix_combination
      (
        storage::Set< util::SiPtr< const SSEGeometryPhiPsi> >( helix_geometries.Begin(), helix_geometries.End()),
        helix_geometries.GetSize() - HELICES
      );
      math::Combination< util::SiPtr< const SSEGeometryPhiPsi> > strand_combination
      (
        storage::Set< util::SiPtr< const SSEGeometryPhiPsi> >( strand_geometries.Begin(), strand_geometries.End()),
        strand_geometries.GetSize() - STRANDS
      );

      // get the possible combinations
      storage::List< storage::Set< util::SiPtr< const SSEGeometryPhiPsi> > > helix_combinations
      (
        helix_combination.GetAllCombinations()
      );
      storage::List< storage::Set< util::SiPtr< const SSEGeometryPhiPsi> > > strand_combinations
      (
        strand_combination.GetAllCombinations()
      );

      // initialize list of potential subdomains
      util::ShPtrList< FoldTemplate> subdomains;

      // iterate through the helix combinations
      for
      (
        storage::List< storage::Set< util::SiPtr< const SSEGeometryPhiPsi> > >::const_iterator
          helix_itr( helix_combinations.Begin()), helix_itr_end( helix_combinations.End());
        helix_itr != helix_itr_end; ++helix_itr
      )
      {
        // iterate through the strand combinations
        for
        (
          storage::List< storage::Set< util::SiPtr< const SSEGeometryPhiPsi> > >::const_iterator
            strand_itr( strand_combinations.Begin()), strand_itr_end( strand_combinations.End());
          strand_itr != strand_itr_end; ++strand_itr
        )
        {
          // copy the helix set
          storage::Set< util::SiPtr< const SSEGeometryPhiPsi> > combined_set( *helix_itr);

          // insert the strand entries
          combined_set.InsertElements( *strand_itr);

          // make a copy of the graph in the topology
          // the topology will have the new vertices but they will still have the old edges
          // this is risky but OK to use here since we remove unused helices and strands to get the sub-topology
          util::ShPtr< Topology::GraphType> sp_graph( m_Topology.GetGraph().HardCopy());

          // initialize geometry vector
          util::ShPtrVector< SSEGeometryPhiPsi> geometry_vector( m_Geometries);

          // iterate over the set
          for
          (
            storage::Set< util::SiPtr< const SSEGeometryPhiPsi> >::const_iterator set_itr( combined_set.Begin()),
              set_itr_end( combined_set.End());
            set_itr != set_itr_end; ++set_itr
          )
          {
            // remove the corresponding vertex
            sp_graph->DeleteVertex( sp_graph->FindVertex( *set_itr));
            geometry_vector.RemoveElements( geometry_vector.Find( *set_itr));
          }

          // initialize bool to check whether any vertices are unconnected
          bool has_unconnected_vertex( false);

          // if the graph has more than one vertex
          if( sp_graph->GetVertices().GetSize() > 1)
          {
            // iterate through the graph vertices
            for
            (
              Topology::GraphType::VertexContainerType::const_iterator vertex_itr( sp_graph->GetVertices().Begin()),
                vertex_itr_end( sp_graph->GetVertices().End());
              vertex_itr != vertex_itr_end; ++vertex_itr
            )
            {
              // if the number of edges is zero, set the bool to true and break
              if( ( *vertex_itr)->GetDegree() == 0)
              {
                has_unconnected_vertex = true;
                break;
              }

              // if the vertex is a strand
              if( STRANDS > 1 && ( *vertex_itr)->GetData()->GetType() == biol::GetSSTypes().STRAND)
              {
                // initialize bool for strand pair
                bool strand_pair( false);

                // iterate over the edges
                for
                (
                  Topology::GraphType::EdgeContainerType::const_iterator edge_itr( ( *vertex_itr)->GetEdges().Begin()),
                    edge_itr_end( ( *vertex_itr)->GetEdges().End());
                  edge_itr != edge_itr_end; ++edge_itr
                )
                {
                  // if the edge is a strand-strand type
                  if( edge_itr->GetData().GetContactType() == contact::GetTypes().STRAND_STRAND)
                  {
                    // found a strand pairing
                    strand_pair = true;
                    break;
                  }
                }

                // if this is a single strand
                if( !strand_pair)
                {
                  // break from iteration
                  has_unconnected_vertex = true;
                  break;
                }
              }
            }
          }

          // if the graph has no unconnected vertices
          if( !has_unconnected_vertex && !sp_graph->GetVertices().IsEmpty())
          {
            subdomains.PushBack
            (
              util::ShPtr< FoldTemplate>( new FoldTemplate( geometry_vector, m_TopologyCollector, m_PDBID))
            );
          }
        }
      }

      BCL_MessageDbg
      (
        "Number of subgraphs: " + util::Format()( subdomains.GetSize()) + " out of " +
        util::Format()( helix_combinations.GetSize() * strand_combinations.GetSize())
      );

      // end
      return subdomains;
    }

  } // namespace assemble
} // namespace bcl
