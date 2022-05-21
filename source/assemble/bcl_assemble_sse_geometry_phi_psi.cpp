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
#include "assemble/bcl_assemble_sse_geometry_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_aa_sequence_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPhiPsi::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPhiPsi())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPhiPsi::SSEGeometryPhiPsi() :
      m_Geometry(),
      m_PhiPsi()
    {
    }

    //! @brief construct from an SSE
    //! @param ORIGINAL_SSE SSE to build from
    SSEGeometryPhiPsi::SSEGeometryPhiPsi( const SSE &ORIGINAL_SSE) :
      m_Geometry( new SSEGeometry( ORIGINAL_SSE)),
      m_PhiPsi( new biol::AASequencePhiPsi( ORIGINAL_SSE))
    {
    }

    //! @brief construct from phi/psi angles
    //! @param PHI_PSI phi/psi angles
    //! @param SS_TYPE SS type to use to build the geometry
    //! @param IDENTIFICATION identification for this geometry
    //! @param SET_GEOMETRY bool whether to set the geometry using the given SS type, if it is coil, or undefined,
    //!        the geometry is not constructed
    SSEGeometryPhiPsi::SSEGeometryPhiPsi
    (
      const biol::AASequencePhiPsi &PHI_PSI,
      const biol::SSType &SS_TYPE,
      const std::string &IDENTIFICATION,
      const bool SET_GEOMETRY
    ) :
      m_Geometry( new SSEGeometry( SS_TYPE, IDENTIFICATION)),
      m_PhiPsi( PHI_PSI.Clone())
    {
      // build the geometry if requested
      if( SET_GEOMETRY && SS_TYPE.IsDefined() && SS_TYPE->IsStructured())
      {
        SetSSEGeometryUsingPhiPsi();
      }
    }

    //! @brief copy constructor
    //! @param GEOMETRY SSEGeometryPhiPsi to be copied
    SSEGeometryPhiPsi::SSEGeometryPhiPsi( const SSEGeometryPhiPsi &GEOMETRY) :
      m_Geometry( GEOMETRY.m_Geometry.HardCopy()),
      m_PhiPsi( GEOMETRY.m_PhiPsi.HardCopy())
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEGeometryPhiPsi
    SSEGeometryPhiPsi *SSEGeometryPhiPsi::Clone() const
    {
      return new SSEGeometryPhiPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEGeometryPhiPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief use the phi/psi information to build the SSEGeometry object
    void SSEGeometryPhiPsi::SetSSEGeometryUsingPhiPsi()
    {
      // construct a poly-A string of a length to match the number of phi/psis
      const std::string dummy_string( m_PhiPsi->GetAngles().GetSize(), 'A');

      // build a dummy sequence using the phi/psi angles
      biol::AASequence dummy_seq
      (
        biol::AASequenceFactory::BuildSequenceFromFASTAString( dummy_string, biol::GetAAClasses().e_AABackBone)
      );

      // fit the dummy sequence
      biol::AASequenceFactory::FitSequence( dummy_seq, *m_PhiPsi, GetType());

      // update the geometry by constructing an SSE and using it to build a SSEGeometry
      const std::string identification( GetIdentification());
      m_Geometry = util::ShPtr< SSEGeometry>( new SSEGeometry( SSE( dummy_seq, GetType())));
      m_Geometry->SetIdentification( identification);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Translation to be applied
    void SSEGeometryPhiPsi::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_Geometry->Translate( TRANSLATION);
      m_PhiPsi->Translate( TRANSLATION);
    }

    //! @brief transform the object by a given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void SSEGeometryPhiPsi::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      m_Geometry->Transform( TRANSFORMATION_MATRIX_3D);
      m_PhiPsi->Transform( TRANSFORMATION_MATRIX_3D);
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void SSEGeometryPhiPsi::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      m_Geometry->Rotate( ROTATION_MATRIX_3D);
      m_PhiPsi->Rotate( ROTATION_MATRIX_3D);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator = SSEGeometryPhiPsi assign operator
    //! @param SSE_GEOMETRY_PHI_PSI_RHS SSEGeometry to be assigned to
    //! @return this SSEGeometryPhiPsi after assignment
    SSEGeometryPhiPsi &SSEGeometryPhiPsi::operator =( const SSEGeometryPhiPsi &SSE_GEOMETRY_PHI_PSI_RHS)
    {
      // assign members
      m_Geometry = SSE_GEOMETRY_PHI_PSI_RHS.m_Geometry.HardCopy();
      m_PhiPsi = SSE_GEOMETRY_PHI_PSI_RHS.m_PhiPsi.HardCopy();

      //end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPhiPsi::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Geometry, ISTREAM);
      io::Serialize::Read( m_PhiPsi, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPhiPsi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Geometry, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PhiPsi, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  /////////////////////////
  // SSEGeometryLessThan //
  /////////////////////////

    //! @brief return true if SSE_GEO_A is less than SSE_GEO_B in size
    //! @param SSE_GEO_A first SSEGeometryPhiPsi
    //! @param SSE_GEO_B second SSEGeometryPhiPsi
    //! @return true if SSE_GEO_A is less than SSE_GEO_B in size
    bool SSEGeometryPhiPsiLessThan::operator()
    (
      const SSEGeometryPhiPsi &SSE_GEO_A,
      const SSEGeometryPhiPsi &SSE_GEO_B
    ) const
    {
      return SSEGeometryInterfaceLessThan()( SSE_GEO_A, SSE_GEO_B);
    }

    //! @brief return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
    //! @param PTR_SSE_GEO_A first SSEGeometryPhiPsi
    //! @param PTR_SSE_GEO_B second SSEGeometryPhiPsi
    //! @return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
    bool SSEGeometryPhiPsiLessThan::operator()
    (
      const util::PtrInterface< SSEGeometryPhiPsi> &PTR_SSE_GEO_A,
      const util::PtrInterface< SSEGeometryPhiPsi> &PTR_SSE_GEO_B
    ) const
    {
      return operator()( *PTR_SSE_GEO_A, *PTR_SSE_GEO_B);
    }

    //! @brief return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
    //! @param PTR_SSE_GEO_A first SSEGeometryPhiPsi
    //! @param PTR_SSE_GEO_B second SSEGeometryPhiPsi
    //! @return true if PTR_SSE_GEO_A is less than PTR_SSE_GEO_B in size
    bool SSEGeometryPhiPsiLessThan::operator()
    (
      const util::PtrInterface< const SSEGeometryPhiPsi> &PTR_SSE_GEO_A,
      const util::PtrInterface< const SSEGeometryPhiPsi> &PTR_SSE_GEO_B
    ) const
    {
      return operator()( *PTR_SSE_GEO_A, *PTR_SSE_GEO_B);
    }

  } // namespace assemble
  
} // namespace bcl
