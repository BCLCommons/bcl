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
#include "score/bcl_score_density_profile_sse_agreement.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_map_cylindrical.h"
#include "restraint/bcl_restraint_body.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  ///////////
  // enums //
  ///////////

    //! @brief Profile1DType as string
    //! @param PROFILE_1D_TYPE the init type
    //! @return the Profile1DType as string
    const std::string &DensityProfileSSEAgreement::GetProfile1DTypeString( const DensityProfileSSEAgreement::Profile1DType &PROFILE_1D_TYPE)
    {
      static const std::string s_profile_1d_type_strings[] =
      {
        "Angle",
        "Height",
        "Radius",
        GetStaticClassName< Profile1DType>()
      };

      return s_profile_1d_type_strings[ PROFILE_1D_TYPE];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DensityProfileSSEAgreement::s_Instance
    (
      GetObjectInstances().AddInstance( new DensityProfileSSEAgreement())
    );

    const size_t DensityProfileSSEAgreement::s_NumberWedges = 36;
    const double DensityProfileSSEAgreement::s_UpperRadius = 10.0;
    const double DensityProfileSSEAgreement::s_LowerRadius = 2.0;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DensityProfileSSEAgreement::DensityProfileSSEAgreement()
    {
    }

    //! @brief construct from Simulator, density map and body restraints
    //! @param SP_SIMULATOR ShPtr to a simulator
    //! @param DENSITY_MAP the experimental density map
    //! @param BODY_RESTRAINT bodies
    //! @param PROFILE_TYPE 1d profile to be used
    DensityProfileSSEAgreement::DensityProfileSSEAgreement
    (
      const util::ShPtr< density::SimulateInterface> &SP_SIMULATOR,
      const density::Map &DENSITY_MAP,
      const restraint::Body &BODY_RESTRAINT,
      const Profile1DType PROFILE_TYPE
    ) :
      m_Profiles(),
      m_Simulator( SP_SIMULATOR),
      m_HeightResolution( DENSITY_MAP.GetCellWidth().X()),
      m_RadiusResolution( DENSITY_MAP.GetCellWidth().Y()),
      m_ProfileType( PROFILE_TYPE)
    {
      // get all bodies/geometries from restraints
      const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &bodies( BODY_RESTRAINT.GetBody());

      // iterate over all geometries
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator body_itr( bodies->Begin()), body_itr_end( bodies->End());
        body_itr != body_itr_end;
        ++body_itr
      )
      {
        // create cylindrical map from current body and given density map
        const density::MapCylindrical density_map
        (
          **body_itr,
          DENSITY_MAP,
          m_HeightResolution,
          m_RadiusResolution,
          s_NumberWedges,
          s_UpperRadius
        );

        // from cylindrical density map, histogram is created representing the profile with option to select histogram type
        math::Histogram profile;
        switch( m_ProfileType)
        {
          case e_Angle:
            profile = density_map.OneDProfileAngle( s_LowerRadius, s_UpperRadius);
            break;
          case e_Height:
            profile = density_map.OneDProfileHeight( s_LowerRadius, s_UpperRadius);
            break;
          case e_Radius:
            profile = density_map.OneDProfileRadius( s_LowerRadius, s_UpperRadius);
            break;
          case s_NumberProfile1DTypes:
            break;
        }

        // insert profile
        m_Profiles[ util::SiPtr< const assemble::SSEGeometryInterface>( *body_itr)] = profile;
      }
    }

    //! @brief Clone function
    //! @return pointer to new DensityProfileSSEAgreement
    DensityProfileSSEAgreement *DensityProfileSSEAgreement::Clone() const
    {
      return new DensityProfileSSEAgreement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DensityProfileSSEAgreement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get scheme
    //! @return scheme
    const std::string &DensityProfileSSEAgreement::GetScheme() const
    {
      return GetProfile1DTypeString( m_ProfileType);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() takes a BODY and an SSE and compares their profiles
    //! @param BODY the density rod that is occupied by the given SSE
    //! @param SECONDAY_STRUCTURE_ELEMENT the secondary structure that was assigned to that density rod
    //! @return return a double which is the agreement of the rod with the SSE
    double DensityProfileSSEAgreement::operator()
    (
      const assemble::SSEGeometryInterface &BODY,
      const assemble::SSE &SECONDAY_STRUCTURE_ELEMENT
    ) const
    {
      // simulate density for that sse using m_Simulator
      const density::Map sse_density_map( m_Simulator->operator()( SECONDAY_STRUCTURE_ELEMENT.GetAtoms()));

      // calculate profile for that sse
      const density::MapCylindrical sse_map_cylindrical
      (
        BODY,
        sse_density_map,
        m_HeightResolution,
        m_RadiusResolution,
        s_NumberWedges,
        s_UpperRadius
      );

      // from SSE, histogram is created representing the profile with option to select histogram type
      math::Histogram sse_profile;
      switch( m_ProfileType)
      {
        case e_Angle:
          sse_profile = sse_map_cylindrical.OneDProfileAngle( s_LowerRadius, s_UpperRadius);
          break;
        case e_Height:
          sse_profile = sse_map_cylindrical.OneDProfileHeight( s_LowerRadius, s_UpperRadius);
          break;
        case e_Radius:
          sse_profile = sse_map_cylindrical.OneDProfileRadius( s_LowerRadius, s_UpperRadius);
          break;
        case s_NumberProfile1DTypes:
          break;
      }

      // find profile for density rod (second or first)
      storage::Map< util::SiPtr< const assemble::SSEGeometryInterface>, math::Histogram>::const_iterator itr( m_Profiles.Find( util::SiPtr< const assemble::SSEGeometryInterface>( BODY)));
      BCL_Assert( itr != m_Profiles.End(), "given density rod was not given to initialization");
      const math::Histogram &rod_profile( itr->second);

      // calculate the ccc
      const double ccc
      (
        CCC
        (
          sse_profile.GetHistogram().Begin(), sse_profile.GetHistogram().End(),
          rod_profile.GetHistogram().Begin(), rod_profile.GetHistogram().End()
        )
      );

      // return the CCC between the profiles
      return -ccc;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DensityProfileSSEAgreement::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Profiles, ISTREAM);
      io::Serialize::Read( m_Simulator, ISTREAM);
      io::Serialize::Read( m_HeightResolution, ISTREAM);
      io::Serialize::Read( m_RadiusResolution, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DensityProfileSSEAgreement::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Profiles        , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Simulator       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HeightResolution, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RadiusResolution, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the CCC between two arrays of doubles
    //! @param BEG_PTR_A
    //! @param END_PTR_A
    //! @param BEG_PTR_B
    //! @param END_PTR_B
    //! @return the cross correlation coefficient between range a and b
    double DensityProfileSSEAgreement::CCC( const double *BEG_PTR_A, const double *END_PTR_A, const double *BEG_PTR_B, const double *END_PTR_B)
    {
      // number of voxels that are above the CONTOUR_LEVEL in the experimental (this) density map
      size_t count_voxel( 0);
      double sum_sim( 0);
      double sum_exp( 0);
      double sum_sim2( 0);
      double sum_exp2( 0);
      double sum_exp_sim( 0);

      // iterate over experimental and simulated sub tensor
      for
      (
        const double *exp( BEG_PTR_A), *exp_end( END_PTR_A), *sim( BEG_PTR_B), *sim_end( END_PTR_B);
        exp != exp_end && sim != sim_end;
        ++exp, ++sim
      )
      {
        const double sim_int( *sim);
        const double exp_int( *exp);

        ++count_voxel;
        sum_exp += exp_int;
        sum_sim += sim_int;
        sum_exp_sim += exp_int * sim_int;
        sum_exp2 += math::Sqr( exp_int);
        sum_sim2 += math::Sqr( sim_int);
      }

      // calculate actual correlation
      double correlation( count_voxel * sum_exp_sim - sum_exp * sum_sim);
      correlation /= math::Sqrt( count_voxel * sum_exp2 - math::Sqr( sum_exp)) * math::Sqrt( count_voxel * sum_sim2 - math::Sqr( sum_sim));

      // end
      return correlation;
    }

  } // namespace score
} // namespace bcl
