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

#ifndef BCL_ASSEMBLE_FOLD_TEMPLATE_H_
#define BCL_ASSEMBLE_FOLD_TEMPLATE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_collector_topology_interface.h"
#include "bcl_assemble_sse_geometry_compare.h"
#include "bcl_assemble_sse_geometry_packing_pickers.h"
#include "bcl_assemble_sse_geometry_phi_psi.h"
#include "bcl_assemble_topology.h"
#include "coord/bcl_coord_movable_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FoldTemplate
    //! @brief This class stores a collection of bodies as a fold template
    //! @details This class describes fold templates, which are collections of bodies describing SSEs
    //!
    //! @see @link example_assemble_fold_template.cpp @endlink
    //! @author weinerbe
    //! @date Dec 9, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FoldTemplate :
      public coord::MovableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Vector of geometries in the fold template
      util::ShPtrVector< SSEGeometryPhiPsi> m_Geometries;

      //! PDB ID of protein fold template was created from
      std::string m_PDBID;

      //! topology representing this fold template
      Topology m_Topology;

      //! collector to use to calculate the topology
      util::ShPtr< CollectorTopologyInterface> m_TopologyCollector;

      //! bool whether the geometries should be sorted by size
      bool m_SortBySize;

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
      FoldTemplate();

      //! @brief construct from a vectors geometries
      //! @param GEOMETRIES bodies of SSEs in template
      //! @param TOPOLOGY_COLLECTOR collector to use to calculate the topology
      //! @param PDB_ID PDB ID of passed protein model
      //! @param SORT_BY_SIZE whether to sort the geometries by size
      FoldTemplate
      (
        const util::ShPtrVector< SSEGeometryPhiPsi> &GEOMETRIES,
        const util::ShPtr< CollectorTopologyInterface> &TOPOLOGY_COLLECTOR,
        const std::string &PDB_ID = GetDefaultPDBID(),
        const bool SORT_BY_SIZE = true
      );

      //! @brief construct from a protein model and PDB ID, sorting the geometries by size
      //! @param PROTEIN_MODEL protein model
      //! @param TOPOLOGY_COLLECTOR collector to use to calculate the topology
      //! @param PDB_ID PDB ID of passed protein model
      //! @param TRIM_GEOMETRIES whether to trim the geometries by 1 residue at each terminus
      FoldTemplate
      (
        const ProteinModel &PROTEIN_MODEL,
        const util::ShPtr< CollectorTopologyInterface> &TOPOLOGY_COLLECTOR,
        const std::string &PDB_ID = GetDefaultPDBID(),
        const bool TRIM_GEOMETRIES = true
      );

      //! @brief copy constructor
      //! @param FOLD_TEMPLATE FoldTemplate to be copied
      FoldTemplate( const FoldTemplate &FOLD_TEMPLATE);

      //! @brief Clone function
      //! @return pointer to new FoldTemplate
      FoldTemplate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief calls GetIndentification on each geometry
      //! @return std::string of GetIndentification from each geometry, separated by a new line
      std::string GetGeometryIndentifications() const;

      //! @brief get the helical geometries
      //! @return the helical geometries
      util::SiPtrVector< const SSEGeometryPhiPsi> GetHelicalGeometries() const;

      //! @brief get the strand geometries
      //! @return the strand geometries
      util::SiPtrVector< const SSEGeometryPhiPsi> GetStrandGeometries() const;

      //! @brief get the geometries that make up the fold template
      //! @return the geometries that make up the fold template
      const util::ShPtrVector< SSEGeometryPhiPsi> &GetGeometries() const
      {
        return m_Geometries;
      }

      //! @brief get the PDB ID of the protein the fold template was created from
      //! @return the PDB ID of the protein the fold template was created from
      const std::string &GetPDBID() const
      {
        return m_PDBID;
      }

      //! @brief returns whether the topology is initialized or not
      //! @return whether the topology is initialized or not
      bool IsTopologyInitialized() const;

      //! @brief creates the topology from the geometry information
      //! @return the topology
      void CalculateTopology();

      //! @brief gets the topology
      //! @return the topology
      const Topology &GetTopology() const
      {
        return m_Topology;
      }

      //! @brief set the topology
      //! @param TOPOLOGY Topology to be copied
      void SetTopology( const Topology &TOPOLOGY);

      //! @brief get the center of mass of the fold template
      //! @return the center of mass of the fold template
      linal::Vector3D GetCenter() const;

      //! @brief get the radius of gyration of the fold template
      //! @return the radius of gyration of the fold template
      double GetRadiusOfGyration() const;

      //! @brief gets a subdomain of the appropriate size
      //! @param HELICES number of helices in the subdomain
      //! @param STRANDS number of strands in the subdomain
      //! @return a subdomain of the appropriate size
      FoldTemplate GetSubDomain( const size_t HELICES, const size_t STRANDS) const;

      //! @brief gets a subdomain that matches the passed SSEs by given comparison
      //! @param SSES sses used to chose subdomains based on length
      //! @param SSE_GEOMETRY_COMPARE comparison method
      //! @return a subdomain that matches the passed SSEs by given comparison
      FoldTemplate GetSubDomain
      (
        const util::SiPtrVector< const SSE> &SSES,
        const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &SSE_GEOMETRY_COMPARE =
          SSEGeometryWithinSizeTolerance()
      ) const;

      //! @brief gets the default PDB ID to be used if none is supplied
      //! @return the default PDB ID to be used if none is supplied
      static const std::string &GetDefaultPDBID();

      //! @brief returns whether the geometries are close in size to the passed SSEs
      //! @param SSES sses to compare geometries to
      //! @param GEOMETRY_COMPARE function to decide whether an SSE should be fit into a geometry
      //! @return whether the geometries are close in size to the passed SSEs
      bool HasSimilarSizeGeometries
      (
        const util::SiPtrVector< const SSE> &SSES,
        const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &GEOMETRY_COMPARE
      ) const;

      //! @brief returns whether there is at least one geometry that is similar to the SSE
      //! @param SP_SSE sse to compare geometries to
      //! @param GEOMETRY_COMPARE function to decide whether an SSE should be fit into a geometry
      //! @return whether there is at least one geometry that is similar to the SSE
      bool HasSimilarSizeGeometry
      (
        const util::SiPtr< const SSE> &SP_SSE,
        const math::BinaryFunctionInterface< SSE, SSEGeometryPhiPsi, bool> &GEOMETRY_COMPARE
      ) const;

      //! @brief returns whether all the geometries are defined
      //! @return whether all the geometries are defined
      bool HasDefinedGeometries() const;

      //! @brief returns whether all the angles are defined
      //! @return whether all the angles are defined
      bool HasDefinedAngles() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief performs translation given by linal::Vector3D
      //! @param TRANSLATION linal::Vector3D to be used for translation
      void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief performs rotation & translation given by 4D transformation matrix
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be used for transformation
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D);

      //! @brief performs rotation given by math::RotationMatrix3D
      //! @param ROTATIONMATRIX3D math::RotationMatrix3D to be used for rotation
      void Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D);

      //! @brief fits the given SSEs into this fold template and returns a domain
      //! @param SSES_TO_FIT vector of SSEs to place into the fold template
      //! @return a domain with the placed SSES
      Domain FitSSEs( const util::SiPtrVector< const SSE> &SSES_TO_FIT) const;

      //! @brief calculates the SSEGeometries using the associated phi/psi information
      void CalculateSSEGeometries();

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      //! @param FOLD_TEMPLATE FoldTemplate to be copied
      //! @return This template after all members are assigned to values from FOLD_TEMPLATE
      FoldTemplate &operator =( const FoldTemplate &FOLD_TEMPLATE);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write minimal information needed to construct a FoldTemplate to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @return outputstream which was written to
      std::ostream &WriteCompact( std::ostream &OSTREAM) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief create geometry member variables
      //! @param PROTEIN_MODEL protein model
      //! @param TRIM_GEOMETRIES whether to trim the geometries
      void CreateGeometries( const ProteinModel &PROTEIN_MODEL, const bool TRIM_GEOMETRIES);

      //! @brief calculate coordinates for each body in a fold template
      //! @return coordinates for a fold template
      storage::Vector< linal::Vector3D> CollectCoordinates() const;

      //! @brief writes the given coordinates to one line
      //! @param COORDS Vector3D to be written
      //! @return string containing the coordinates on one line
      static std::string WriteCoordinates( const linal::Vector3D &COORDS);

      //! @brief trims vectors of SSEs and geometries (whichever is biggest so that they match)
      //! @param SSES vector of SSEs
      //! @param GEOS vector of SSEGeometries
      static void MatchSSEsAndGeometries
      (
        util::ShPtrVector< SSE> &SSES, util::SiPtrVector< const SSEGeometryPhiPsi> &GEOS
      );

      //! @brief get the average number of close geometries
      //! @return the average number of close geometries
      double AverageNumberOfCloseGeometries() const;

      //! @brief get the number of geometries within DISTANCE
      //! @param GEOMETRY geometry of interest
      //! @return the number of bodies within DISTANCE
      size_t NumberOfCloseGeometries( const util::ShPtr< SSEGeometryPhiPsi> &GEOMETRY) const;

      //! @brief get geometries that contact the passed geometry in the fold template
      //! @param GEOMETRY geometry of interest
      //! @return ShPtrVector of geometries that contact the passed geometry in the fold template
      util::ShPtrVector< SSEGeometryPhiPsi> GetNeighborGeometries( const util::ShPtr< SSEGeometryPhiPsi> &GEOMETRY) const;

      //! @brief gets all the possible subdomains of the requested size
      //! @param HELICES number of helices in the subdomain
      //! @param STRANDS number of strands in the subdomain
      //! @return all the possible subdomains of the requested size
      util::ShPtrList< FoldTemplate> GetAllSubdomains( const size_t HELICES, const size_t STRANDS) const;

    }; // class FoldTemplate

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_FOLD_TEMPLATE_H_ 
