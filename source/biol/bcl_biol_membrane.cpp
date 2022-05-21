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
#include "biol/bcl_biol_membrane.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_environment_types.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! @brief return command line flag for defining the membrane thickness
    //! @return command line flag for defining the membrane thickness
    util::ShPtr< command::FlagInterface> &Membrane::GetFlagMembrane()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic( "membrane", "flag for using a membrane and defining membrane thicknesses")
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( GetParameterCoreThickness());
        flag->PushBack( GetParameterTransitionThickness());
        flag->PushBack( GetParameterGapThickness());
      }

      // end
      return s_flag;
    }

    //! @brief return command line parameter for defining the membrane core thickness
    //! @return command line parameter for defining the membrane core thickness
    util::ShPtr< command::ParameterInterface> &Membrane::GetParameterCoreThickness()
    {
      // initialize static instance of the parameter
      static util::ShPtr< command::ParameterInterface> s_core_thickness_param
      (
        new command::Parameter
        (
          "membrane_core",
          "half the thickness of membrane core region in Angstroem",
          command::ParameterCheckRanged< double>( 0.0, 50.0),
          "10.0"
        )
      );

      // end
      return s_core_thickness_param;
    }

    //! @brief return command line parameter for defining the membrane transition thickness
    //! @return command line parameter for defining the membrane transition thickness
    util::ShPtr< command::ParameterInterface> &Membrane::GetParameterTransitionThickness()
    {
      // initialize static instance of the parameter
      static util::ShPtr< command::ParameterInterface> s_transition_thickness_param
      (
        new command::Parameter
        (
          "membrane_trans",
          "thickness of membrane transition region in Angstroem",
          command::ParameterCheckRanged< double>( 0.0, 50.0),
          "10.0"
        )
      );

      // end
      return s_transition_thickness_param;
    }

    //! @brief return command line parameter for defining the membrane gap thickness
    //! @return command line parameter for defining the membrane gap thickness
    util::ShPtr< command::ParameterInterface> &Membrane::GetParameterGapThickness()
    {
      // initialize static instance of the parameter
      static util::ShPtr< command::ParameterInterface> s_gap_thickness_param
      (
        new command::Parameter
        (
          "membrane_gap",
          "thickness of membrane gap between regions in Angstroem",
          command::ParameterCheckRanged< double>( 0.0, 50.0),
          "2.5"
        )
      );

      // end
      return s_gap_thickness_param;
    }

    //! @brief create a membrane object from commandline arguments
    //! @return membrane object created from commandline arguments
    Membrane Membrane::GetCommandLineMembrane()
    {
      // construct and return a membrane object constructed from commandline arguments
      return Membrane
      (
        GetParameterCoreThickness()->GetNumericalValue< double>(),
        GetParameterTransitionThickness()->GetNumericalValue< double>(),
        GetParameterGapThickness()->GetNumericalValue< double>()
      );
    }

    //! @brief static undefined membrane object
    //! @return undefined membrane object
    const Membrane &Membrane::GetUndefinedMembrane()
    {
      static const Membrane s_membrane
      (
        util::GetUndefined< double>(),
        util::GetUndefined< double>(),
        util::GetUndefined< double>()
      );
      return s_membrane;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Membrane::s_Instance
    (
      GetObjectInstances().AddInstance( new Membrane())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Membrane::Membrane() :
      m_Orientation(),
      m_Thicknesses
      (
        FillThicknessVector
        (
          GetEnvironmentTypes().e_MembraneCore->GetDefaultThickness(),
          GetEnvironmentTypes().e_Transition->GetDefaultThickness(),
          GetEnvironmentTypes().e_GapCoreTransition->GetDefaultThickness()
        )
      ),
      m_Limits( FillLimitsVector( m_Thicknesses))
    {
      m_IsDefined = m_Orientation.IsDefined() && IsDefined( m_Thicknesses) && IsDefined( m_Limits);
    }

    //! @brief constructor from membrane normal, all thicknesses and gap thickness
    //! @param THICKNESSES is a vector of membrane thicknesses
    //! @param CENTER center of the membrane
    //! @param NORMAL membrane normal
    Membrane::Membrane
    (
      const storage::Vector< double> &THICKNESSES,
      const linal::Vector3D &NORMAL,
      const linal::Vector3D &CENTER
    ) :
      m_Orientation( OrientationFromNormal( NORMAL, CENTER)),
      m_Thicknesses( THICKNESSES),
      m_Limits( FillLimitsVector( m_Thicknesses))
    {
      m_IsDefined = m_Orientation.IsDefined() && IsDefined( m_Thicknesses) && IsDefined( m_Limits);
      BCL_Assert
      (
        THICKNESSES.GetSize() <= GetEnvironmentTypes().GetEnumCount(),
        "given vector of thicknesses had too many elements"
      );
    }

    //! @brief constructor from membrane normal, specified thicknesses and gap thickness
    //! @param THICKNESS_CORE thickness of membrane core
    //! @param THICKNESS_TRANSITION thickness of membrane transition region
    //! @param THICKNESS_GAP thickness of membrane gap
    //! @param NORMAL Membrane normal
    Membrane::Membrane
    (
      const double THICKNESS_CORE,
      const double THICKNESS_TRANSITION,
      const double THICKNESS_GAP,
      const linal::Vector3D &NORMAL
    ) :
      m_Orientation( OrientationFromNormal( NORMAL)),
      m_Thicknesses( FillThicknessVector( THICKNESS_CORE, THICKNESS_TRANSITION, THICKNESS_GAP)),
      m_Limits( FillLimitsVector( m_Thicknesses))
    {
      m_IsDefined = m_Orientation.IsDefined() && IsDefined( m_Thicknesses) && IsDefined( m_Limits);
    }

    //! @brief virtual copy constructor
    Membrane *Membrane::Clone() const
    {
      return new Membrane( *this);
    }

    //! @brief destructor
    Membrane::~Membrane()
    {
      m_DestructorSignal.Emit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Membrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the membrane normal
    //! @return the membrane normal
    linal::Vector3D Membrane::GetNormal() const
    {
      linal::Vector3D normal( GetAxis( coord::GetAxes().e_Z));
      normal.Normalize();
      return normal;
    }

    //! @brief return the Thickness of the region
    //! @return the Thickness of the region
    double Membrane::GetThickness( const EnvironmentType &ENVIRONMENT) const
    {
      return m_Thicknesses( ENVIRONMENT);
    }

    //! @brief return the Limit of the region
    //! @return the Limit of the region
    double Membrane::GetLimit( const EnvironmentType &ENVIRONMENT) const
    {
      return m_Limits( ENVIRONMENT);
    }

    //! @brief returns true if the membrane is not set to static undefined membrane
    //! @return true if the membrane is not set to static undefined membrane
    bool Membrane::IsDefined() const
    {
      return m_IsDefined;
    }

    //! @brief returns the geometric center of the object
    //! @return the geometric center of the object
    linal::Vector3D Membrane::GetCenter() const
    {
      return m_Orientation.GetOrigin();
    }

    //! @brief return the orientation of the object
    //! @return orientation
    linal::Vector3D Membrane::GetAxis( const coord::Axis &AXIS) const
    {
      return m_Orientation.GetAxis( AXIS);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns EnvironmentType according to the coordinates
    //! @param COORDINATES coordinates of object
    //! @return EnvironmentType according to the coordinates
    EnvironmentType Membrane::DetermineEnvironmentType( const linal::Vector3D &COORDINATES) const
    {
      // return undefined environment for undefined coords
      if( !COORDINATES.IsDefined())
      {
        return GetEnvironmentTypes().e_Undefined;
      }

      // calculate the distance from the membrane plane
      const double distance( DistanceFromPlane( COORDINATES));

      //check in which region the distance is the first time smaller than the current limit
      for
      (
        EnvironmentTypes::const_iterator itr( GetEnvironmentTypes().Begin()), itr_end( GetEnvironmentTypes().End());
        itr != itr_end; ++itr
      )
      {
        if( distance < m_Limits( *itr))
        {
          return *itr;
        }
      }

      //if it does not fit in any region, return undefined
      return util::GetUndefined< EnvironmentType>();
    }

    //! @brief returns EnvironmentType according to the z-coordiante and weight, if it is in a gap reagion
    //! @param COORDINATES coordinates of object
    //! @return pair of EnvironmentType and a weight between 0 and one, 1 if it is closer to the innermose region, 0 if it is close to the outer most region
    storage::Pair< EnvironmentType, double>
    Membrane::DetermineEnvironmentTypeAndWeight( const linal::Vector3D &COORDINATES) const
    {
      //determine the absolute z-coordinate
      const double distance( DistanceFromPlane( COORDINATES));

      storage::Pair< EnvironmentType, double> environment_weight( DetermineEnvironmentType( COORDINATES), 1.0);
      const EnvironmentType env_type( environment_weight.First());

      // switch over environment types and change the weight, if a gap gets hit
      if( env_type == GetEnvironmentTypes().e_GapCoreTransition)
      {
        const double position_in_gap
        (
          ( distance - m_Limits( GetEnvironmentTypes().e_MembraneCore)) /
            m_Thicknesses( GetEnvironmentTypes().e_GapCoreTransition) * math::g_Pi
        );
        environment_weight.Second() = math::WeightBetweenZeroAndPi( position_in_gap);
      }
      else if( env_type == GetEnvironmentTypes().e_GapTransitionSolution)
      {
        const double position_in_gap
        (
          ( distance - m_Limits( GetEnvironmentTypes().e_Transition)) /
            m_Thicknesses( GetEnvironmentTypes().e_GapTransitionSolution) * math::g_Pi);
        environment_weight.Second() = math::WeightBetweenZeroAndPi( position_in_gap);
      }
      // if undefined
      else if( !env_type.IsDefined())
      {
        environment_weight.Second() = util::GetUndefined< double>();
      }

      // return determined environment and weight
      return environment_weight;
    }

    //! @brief  returns solvation energy for given z-coordinate and a vector containing three state solvation energies
    //! @param COORDINATES coordinates of object
    //! @return solvation energy for given z-coordinate and a vector containing three state solvation energies
    double Membrane::CalculateSolvationEnergy( const linal::Vector3D &COORDINATES, const linal::Vector3D &TFE) const
    {
      // determine environment and weight
      const storage::Pair< EnvironmentType, double> env_weight( DetermineEnvironmentTypeAndWeight( COORDINATES));

      // if env type is undefined
      if( !env_weight.First().IsDefined())
      {
        return util::GetUndefined< double>();
      }

      // for actual region
      if( !env_weight.First()->IsGap())
      {
        return TFE( env_weight.First()->GetReducedIndex());
      }

      // for gap
      return           env_weight.Second()  * TFE( env_weight.First()->GetReducedIndex() - 1)
             + ( 1.0 - env_weight.Second()) * TFE( env_weight.First()->GetReducedIndex());
    }

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Translation to be applied
    void Membrane::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_Orientation( TRANSLATION);
    }

    //! @brief transform the object by a given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void Membrane::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      m_Orientation( TRANSFORMATION_MATRIX_3D);
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void Membrane::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      m_Orientation( ROTATION_MATRIX_3D);
    }

    //! @brief fills the Thicknesses with given values
    //! @param THICKNESS_CORE thickness of membrane core region
    //! @param THICKNESS_TRANSITION thickness of membrane transition region
    //! @param THICKNESS_GAP thickness of membrane gap region
    //! @return Thickness vector
    storage::Vector< double> Membrane::FillThicknessVector
    (
      const double THICKNESS_CORE,
      const double THICKNESS_TRANSITION,
      const double THICKNESS_GAP
    )
    {
      storage::Vector< double> thicknesses( GetEnvironmentTypes().GetEnumCount());
      thicknesses( GetEnvironmentTypes().e_MembraneCore)          = THICKNESS_CORE;
      thicknesses( GetEnvironmentTypes().e_GapCoreTransition)     = THICKNESS_GAP;
      thicknesses( GetEnvironmentTypes().e_Transition)            = THICKNESS_TRANSITION;
      thicknesses( GetEnvironmentTypes().e_GapTransitionSolution) = THICKNESS_GAP;
      thicknesses( GetEnvironmentTypes().e_Solution)              = GetEnvironmentTypes().e_Solution->GetDefaultThickness();
      thicknesses( GetEnvironmentTypes().e_SolutionInside)        = 0.0;
      thicknesses( GetEnvironmentTypes().e_SolutionOutside)       = 0.0;

      return thicknesses;
    }

    //! @brief Calculates the Limits for the membrane regions
    //! @param THICKNESSES vector of thickness for membrane regions
    //! @return limits for the membrane regions
    storage::Vector< double> Membrane::FillLimitsVector( const storage::Vector< double> &THICKNESSES)
    {
      // vector of limits
      storage::Vector< double> limits( GetEnvironmentTypes().GetEnumCount());

      // set the value
      double sum_thickness( 0.0);

      // for all other limits the limit is the previous limit + the current thickness
      for
      (
        EnvironmentTypes::const_iterator itr( GetEnvironmentTypes().Begin()),
          itr_end( GetEnvironmentTypes().End());
        itr != itr_end; ++itr
      )
      {
        if( itr->GetIndex() < THICKNESSES.GetSize())
        {
          sum_thickness += THICKNESSES( *itr);
        }
        limits( *itr) = sum_thickness;
      }

      return limits;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Membrane::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Orientation, ISTREAM);
      io::Serialize::Read( m_Thicknesses, ISTREAM);
      io::Serialize::Read( m_Limits, ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &Membrane::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write member
      io::Serialize::Write( m_Orientation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Thicknesses, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Limits, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    math::TransformationMatrix3D TransformationMatrixFromPDBTMXML_TMATRIX( std::istream &ISTREAM)
    {
      // create linebuffer
      std::string line_buffer;

      linal::Matrix< double> transformation( 4, 4, double( 0.0));

      std::getline( ISTREAM, line_buffer);
      storage::Vector< std::string> xyz( util::SplitString( line_buffer, "="));
      BCL_Assert( util::TrimString( xyz( 0)) == "<ROWX X", "improper format");
      transformation( 0, 0) = util::ConvertStringToNumericalValue< double>( xyz( 1).substr( 1, xyz( 1).length() - 4));
      transformation( 1, 0) = util::ConvertStringToNumericalValue< double>( xyz( 2).substr( 1, xyz( 2).length() - 4));
      transformation( 2, 0) = util::ConvertStringToNumericalValue< double>( xyz( 3).substr( 1, xyz( 3).length() - 4));
      transformation( 3, 0) = util::ConvertStringToNumericalValue< double>( xyz( 4).substr( 1, xyz( 4).length() - 4));
      std::getline( ISTREAM, line_buffer);
      xyz = util::SplitString( line_buffer, "=");
      BCL_Assert( util::TrimString( xyz( 0)) == "<ROWY X", "improper format");
      transformation( 0, 1) = util::ConvertStringToNumericalValue< double>( xyz( 1).substr( 1, xyz( 1).length() - 4));
      transformation( 1, 1) = util::ConvertStringToNumericalValue< double>( xyz( 2).substr( 1, xyz( 2).length() - 4));
      transformation( 2, 1) = util::ConvertStringToNumericalValue< double>( xyz( 3).substr( 1, xyz( 3).length() - 4));
      transformation( 3, 1) = util::ConvertStringToNumericalValue< double>( xyz( 4).substr( 1, xyz( 4).length() - 4));
      std::getline( ISTREAM, line_buffer);
      xyz = util::SplitString( line_buffer, "=");
      BCL_Assert( util::TrimString( xyz( 0)) == "<ROWZ X", "improper format");
      transformation( 0, 2) = util::ConvertStringToNumericalValue< double>( xyz( 1).substr( 1, xyz( 1).length() - 4));
      transformation( 1, 2) = util::ConvertStringToNumericalValue< double>( xyz( 2).substr( 1, xyz( 2).length() - 4));
      transformation( 2, 2) = util::ConvertStringToNumericalValue< double>( xyz( 3).substr( 1, xyz( 3).length() - 4));
      transformation( 3, 2) = util::ConvertStringToNumericalValue< double>( xyz( 4).substr( 1, xyz( 4).length() - 4));

      // check for valid end
      std::getline( ISTREAM, line_buffer);
      BCL_Assert( util::TrimString( line_buffer) == "</TMATRIX>", "improper format");

      // end
      return math::TransformationMatrix3D( transformation);
    }

    //! @brief create membrane object and transformation matrix from given pdbtm xml file
    //! only the core thickness can be retrieved from the xml file, so that transition region and gap are passed
    //! @param ISTREAM input stream
    //! @param THICKNESS_TRANSITION thickness of the membrane transition region
    //! @param THICKNESS_GAP thickness of the gaps between the regions
    //! @return pair of membrane and transformation matrix
    storage::Pair< Membrane, math::TransformationMatrix3D> Membrane::MembraneAndTransformationFromPDBTMXML
    (
      std::istream &ISTREAM,
      const double THICKNESS_TRANSITION,
      const double THICKNESS_GAP
    )
    {
      // create membrane and transformation matrix, initialize with undefined values
      storage::Pair< Membrane, math::TransformationMatrix3D> membrane_matrix
      (
        Membrane( storage::Vector< double>( GetEnvironmentTypes().GetEnumCount(), util::GetUndefined< double>()), linal::Vector3D()),
        math::TransformationMatrix3D( util::UndefinedObject())
      );

      // create linebuffer
      std::string line_buffer;
      while( std::getline( ISTREAM, line_buffer))
      {
        if( util::TrimString( line_buffer) == "<MEMBRANE>")
        {
          break;
        }
      }

      // if no membrane tag was found, then presumably the protein was soluble, return an undefined membrane
      if( !ISTREAM.good())
      {
        return membrane_matrix;
      }

      // get membrane normale
      std::getline( ISTREAM, line_buffer);
      storage::Vector< std::string> xyz( util::SplitString( line_buffer, "="));
      linal::Vector3D normale;
      BCL_Assert( util::TrimString( xyz( 0)) == "<NORMAL X", "improper format, expected \"<NORMAL X\" but found " + util::TrimString( xyz( 0)));
      normale.X() = util::ConvertStringToNumericalValue< double>( xyz( 1).substr( 1, xyz( 1).length() - 4));
      normale.Y() = util::ConvertStringToNumericalValue< double>( xyz( 2).substr( 1, xyz( 2).length() - 4));
      normale.Z() = util::ConvertStringToNumericalValue< double>( xyz( 3).substr( 1, xyz( 3).length() - 4));

      // get and generate matrix
      std::getline( ISTREAM, line_buffer);
      BCL_Assert( util::TrimString( line_buffer) == "<TMATRIX>", "improper format");

      // prepare
      membrane_matrix.First() = Membrane( normale.Norm(), THICKNESS_TRANSITION, THICKNESS_GAP, normale);
      membrane_matrix.Second() = TransformationMatrixFromPDBTMXML_TMATRIX( ISTREAM);

      // check for valid end
      std::getline( ISTREAM, line_buffer);
      BCL_Assert( util::TrimString( line_buffer) == "</MEMBRANE>", "improper format");

      // some debug information
      BCL_MessageDbg( "transformation matrix found in pdbtm xml file:\n" + util::Format()( membrane_matrix.Second()));

      // end
      return membrane_matrix;
    }

    //! @brief returns pairs of chain id and transformation matrix that has to be applied to the chain to get one monomer for
    //! the biological relevant unit
    //! @param ISTREAM input stream of pdbtm xml file
    //! @param CHAIN_IDS chain ids in the original pdb
    //! @return vector of chain id the transformation is applied to, the new chain id and the matrix that needs to be applied order to get the relevant BIOMOLECULE
    storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >
    Membrane::BioTransformationMatricesFromPDBTMXML( std::istream &ISTREAM, const std::string &CHAIN_IDS)
    {
      std::string delete_chain_ids;
      storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > bio_transformations;

      // create linebuffer and progress to BIOMATRIX
      std::string line_buffer;
      while( std::getline( ISTREAM, line_buffer) && util::TrimString( line_buffer) != "<BIOMATRIX>");

      while( util::TrimString( line_buffer) != "</BIOMATRIX>")
      {
        // error in getting next line
        if( !std::getline( ISTREAM, line_buffer).good())
        {
          break;
        }
        // split line
        storage::Vector< std::string> split_identifier( util::SplitString( line_buffer, "="));

        // check if chain delete was given
        if( util::TrimString( split_identifier( 0)) == "<DELETE CHAINID")
        {
          delete_chain_ids.push_back( split_identifier( 1)[ 1]);
          continue;
        }

        // check if their is a new matrix
        if( util::TrimString( split_identifier( 0)) == "<MATRIX ID")
        {
          // gather all chain ids that this needs to be applied to with its new chain id
          storage::Vector< storage::Pair< char, char> > chain_ids;

          std::getline( ISTREAM, line_buffer);
          // split line
          storage::Vector< std::string> split_apply( util::SplitString( line_buffer, "="));
          while( util::TrimString( split_apply( 0)) == "<APPLY_TO_CHAIN CHAINID")
          {
            const storage::Pair< char, char> chain_id_pair( split_apply( 1)[ 1], split_apply( 2)[ 1]);

            // if the chain id is present in the model
            if( CHAIN_IDS.find( chain_id_pair.First()) != std::string::npos)
            {
              // add the pair
              chain_ids.PushBack( chain_id_pair);
            }
            std::getline( ISTREAM, line_buffer);
            split_apply = util::SplitString( line_buffer, "=");
          }

          // get the matrix for those chains
          BCL_Assert
          (
            util::TrimString( split_apply( 0)) == "<TMATRIX>",
            "improper format. expected <TMATRIX> found: " + split_apply( 0)
          );
          const math::TransformationMatrix3D current_transformation( TransformationMatrixFromPDBTMXML_TMATRIX( ISTREAM));
          std::getline( ISTREAM, line_buffer);
          BCL_Assert
          (
            util::TrimString( line_buffer) == "</MATRIX>",
            "improper format. expected </MATRIX> found: " + line_buffer
          );

          for( storage::Vector< storage::Pair< char, char> >::const_iterator itr( chain_ids.Begin()), itr_end( chain_ids.End()); itr != itr_end; ++itr)
          {
            bio_transformations.PushBack( storage::Triplet< char, char, math::TransformationMatrix3D>( itr->First(), itr->Second(), current_transformation));
          }
        }
      }

      // no bio transformations were found and no chain to delete
      if( bio_transformations.IsEmpty() && delete_chain_ids.empty())
      {
        return bio_transformations;
      }

      // add identity to all given chains, since they are not given in the pdbtmxml file
      for( std::string::const_iterator itr( CHAIN_IDS.begin()), itr_end( CHAIN_IDS.end()); itr != itr_end; ++itr)
      {
        // identify matrix; old chain id is new chain id
        bio_transformations.PushBack( storage::Triplet< char, char, math::TransformationMatrix3D>( *itr, *itr, math::TransformationMatrix3D()));
      }

      // delete undesired chains
      for( std::string::const_iterator char_itr( delete_chain_ids.begin()), char_itr_end( delete_chain_ids.end()); char_itr != char_itr_end; ++char_itr)
      {
        for
        (
          storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> >::iterator
            itr( bio_transformations.Begin()), itr_end( bio_transformations.End()); itr != itr_end; ++itr
        )
        {
          if( itr->Second() == *char_itr)
          {
            bio_transformations.Remove( itr);
            break;
          }
        }
      }

      // end
      return bio_transformations;
    }

    //! @brief evaluates whether all entries in the vector are defined
    //! @param VECTOR vector to be evaulated
    //! @return whether all entries in the vector are defined
    bool Membrane::IsDefined( const storage::Vector< double> &VECTOR)
    {
      // iterate through the vector
      for
      (
        storage::Vector< double>::const_iterator itr( VECTOR.Begin()), itr_end( VECTOR.End());
        itr != itr_end; ++itr
      )
      {
        // if the value is not defined
        if( !util::IsDefined( *itr))
        {
          return false;
        }
      }

      // if this point is reached all values are defined
      return true;
    }

    //! @brief calculates the distance of a point from the center of the membrane plane
    //! @param COORDINATES coordinates to measure
    //! @return calculated distance
    double Membrane::DistanceFromPlane( const linal::Vector3D &COORDINATES) const
    {
      return math::Absolute
      (
        linal::ScalarProduct( GetNormal(), COORDINATES - GetCenter())
      );
    }

    //! @brief gets the membrane orientation from the given normal
    //! @param NORMAL membrane normal
    //! @param CENTER membrane center
    //! @return membrane orientation
    math::TransformationMatrix3D Membrane::OrientationFromNormal
    (
      const linal::Vector3D &NORMAL,
      const linal::Vector3D &CENTER
    )
    {
      // initialize transformation
      math::TransformationMatrix3D transform;

      // get the axis and angle
      linal::Vector3D normal( NORMAL);
      normal.Normalize();
      const linal::Vector3D axis( linal::CrossProduct( normal, coord::GetAxes().e_Z));
      const double angle( ProjAngle( normal, coord::GetAxes().e_Z));

      // only change the orientation if the axis and angle are not zero
      if( axis.Norm() != 0.0 && angle != 0.0)
      {
        transform( math::RotationMatrix3D( axis, angle));
      }

      // apply translation
      transform( CENTER);

      // end
      return transform;
    }

  } // namespace biol
} // namespace bcl
