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
#include "opencl/bcl_opencl_saxs_debye.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math.h"
#include "math/bcl_math_sum_function.h"
#include "opencl/bcl_opencl_vector.h"
#include "restraint/bcl_restraint_sas_data_parameters.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_data.h"
#include "restraint/bcl_restraint_sas_scattering_data.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_stopwatch.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include <CL/cl_platform.h>

namespace bcl
{
  namespace opencl
  {

    class BCL_API FormFactorParameters
    {

    private:

    //////////
    // data //
    //////////

      //! float (c1 parameter)
      float m_ExcludedVolumeParameter;

      //! float (c2 parameter)
      float m_HydrationShellParameter;

      //! float ( sasa)
      float m_SolventAccessableSurfaceArea;

      //! float DisplacedSolventVolume
      float m_DisplacedSolventVolume;

      //! float VanderWaals Radius
      float m_Radius;

      //! float Bound Hydrogen
      float m_BoundHydrogen;

      //! float Crommer Mann Coefficient A1
      float m_A1;

      //! float Crommer Mann Coefficient A2
      float m_A2;

      //! float Crommer Mann Coefficient A3
      float m_A3;

      //! float Crommer Mann Coefficient A4
      float m_A4;

      //! float Crommer Mann Coefficient B1
      float m_B1;

      //! float Crommer Mann Coefficient B2
      float m_B2;

      //! float Crommer Mann Coefficient B3
      float m_B3;

      //! float Crommer Mann Coefficient B4
      float m_B4;

      //! float Crommer Mann Coefficient C
      float m_C;

      //! float X Coordinate
      float m_X;

      //! float Y Coordinate
      float m_Y;

      //! float Z Coordinate
      float m_Z;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @breif Default Constructor
      FormFactorParameters();

      //! @brief Clone function
      //! @return pointer to new FormFactorParameters
      FormFactorParameters *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief set ExcludedVolumeParameter (c1) member variable
      void SetExcludedVolumeParameter( const float &C1)
      {
        m_ExcludedVolumeParameter = C1;
      }

      //! @brief set HydrationShellParameter (c2) member variable
      void SetHydrationShellParameter( const float &C2)
      {
        m_HydrationShellParameter = C2;
      }

      //! @brief set SolventAccessableSurfaceArea ( sasa) member variable
      void SetSolventAccessableSurfaceArea( const double &SASA)
      {
        m_SolventAccessableSurfaceArea = ( float)SASA;
      }

      //! @brief set DisplacedSolventVolume member variable
      void SetDisplacedSolventVolume( const double &DSV)
      {
        m_DisplacedSolventVolume = ( float)DSV;
      }

      //! @brief set VanderWaals Radius member variable
      void SetRadius( const double &RADIUS)
      {
        m_Radius = ( float)RADIUS;
      }

      //! @brief set Bound Hydrogen member variable
      void SetBoundHydrogen( const size_t &BOUND_HYDROGEN)
      {
        m_BoundHydrogen = ( float)BOUND_HYDROGEN;
      }

      //! @brief set Crommer Mann Coefficient A1 member variable
      void SetA1( const double &A1)
      {
        m_A1 = ( float)A1;
      }

      //! @brief set Crommer Mann Coefficient A2 member variable
      void SetA2( const double &A2)
      {
        m_A2 = ( float)A2;
      }

      //! @brief set Crommer Mann Coefficient A3 member variable
      void SetA3( const double &A3)
      {
        m_A3 = ( float)A3;
      }

      //! @brief set Crommer Mann Coefficient A4 member variable
      void SetA4( const double &A4)
      {
        m_A4 = ( float)A4;
      }

      //! @brief set Crommer Mann Coefficient B1 member variable
      void SetB1( const double &B1)
      {
        m_B1 = ( float)B1;
      }

      //! @brief set Crommer Mann Coefficient B2 member variable
      void SetB2( const double &B2)
      {
        m_B2 = ( float)B2;
      }

      //! @brief set Crommer Mann Coefficient B3 member variable
      void SetB3( const double &B3)
      {
        m_B3 = ( float)B3;
      }

      //! @brief set Crommer Mann Coefficient B4 member variable
      void SetB4( const double &B4)
      {
        m_B4 = ( float)B4;
      }

      //! @brief set Crommer Mann Coefficient C member variable
      void SetC( const double &C)
      {
        m_C = ( float)C;
      }

      //! @brief set X Coordinate
      void SetX( const double &X)
      {
        m_X = ( float)X;
      }

       //! @brief set Y Coordinate
      void SetY( const double &Y)
      {
        m_Y = ( float)Y;
      }

      //! @brief set X Coordinate
      void SetZ( const double &Z)
      {
        m_Z = ( float)Z;
      }

      void SetParameters
      (
        const storage::Vector< std::string> &ATOMGROUP,
        const storage::Vector< linal::Vector3D> &COORDINATES,
        const storage::Vector< double> &SASA_VALUE,
        size_t LOCATION,
        float EXCLUDED_VOLUME,
        float HYDRATION_SHELL
      );

      void ShowValues();

      cl_float16 GetParametersAsFloat16()
      {
        cl_float16 host;

        host.s[0] = m_ExcludedVolumeParameter;
        host.s[1] = m_HydrationShellParameter;
        host.s[2] = m_SolventAccessableSurfaceArea;
        host.s[3] = m_DisplacedSolventVolume;
        host.s[4] = m_Radius;
        host.s[5] = m_BoundHydrogen;
        host.s[6] = m_A1;
        host.s[7] = m_A2;
        host.s[8] = m_A3;
        host.s[9] = m_A4;
        host.s[10] = m_B1;
        host.s[11] = m_B2;
        host.s[12] = m_B3;
        host.s[13] = m_B4;
        host.s[14] = m_C;
        host.s[15] = float( 0.0);

        return host;

      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class FormFactorParameters

    FormFactorParameters::FormFactorParameters() :
      m_ExcludedVolumeParameter( util::GetUndefined< float>()),
      m_HydrationShellParameter( util::GetUndefined< float>()),
      m_SolventAccessableSurfaceArea( util::GetUndefined< float>()),
      m_DisplacedSolventVolume( util::GetUndefined< float>()),
      m_Radius( util::GetUndefined< float>()),
      m_BoundHydrogen( util::GetUndefined< float>()),
      m_A1( util::GetUndefined< float>()),
      m_A2( util::GetUndefined< float>()),
      m_A3( util::GetUndefined< float>()),
      m_A4( util::GetUndefined< float>()),
      m_B1( util::GetUndefined< float>()),
      m_B2( util::GetUndefined< float>()),
      m_B3( util::GetUndefined< float>()),
      m_B4( util::GetUndefined< float>()),
      m_C( util::GetUndefined< float>()),
      m_X( util::GetUndefined< float>()),
      m_Y( util::GetUndefined< float>()),
      m_Z( util::GetUndefined< float>())
    {
    }

    void FormFactorParameters::SetParameters
    (
      const storage::Vector< std::string> &ATOMGROUP,
      const storage::Vector< linal::Vector3D> &COORDINATES,
      const storage::Vector< double> &SASA_VALUE,
      size_t LOCATION,
      float EXCLUDED_VOLUME,
      float HYDRATION_SHELL
    )
    {
      SetExcludedVolumeParameter( EXCLUDED_VOLUME);
      SetHydrationShellParameter( HYDRATION_SHELL);
      SetSolventAccessableSurfaceArea( SASA_VALUE( LOCATION));
      SetX( COORDINATES( LOCATION).X());
      SetY( COORDINATES( LOCATION).Y());
      SetZ( COORDINATES( LOCATION).Z());

      if( ATOMGROUP( LOCATION) == "H")
      {
        SetDisplacedSolventVolume( 5.15);
        SetRadius( 1.07);
        SetBoundHydrogen( 0);
        SetA1( 0.493002);
        SetA2( 0.322912);
        SetA3( 0.140191);
        SetA4( 0.040810);
        SetB1( 10.510900);
        SetB2( 26.125700);
        SetB3(  3.142360);
        SetB4( 57.799700);
        SetC( 0.003038);
      }
      else if
      (
          ATOMGROUP( LOCATION) == "C" ||
          ATOMGROUP( LOCATION) == "CH" ||
          ATOMGROUP( LOCATION) == "CH2" ||
          ATOMGROUP( LOCATION) == "CH3"
      )
      {
        SetA1( 2.310000);
        SetA2( 1.020000);
        SetA3( 1.588600);
        SetA4( 0.865000);
        SetB1( 20.843900);
        SetB2( 10.207500);
        SetB3( 0.568700);
        SetB4( 51.651200);
        SetC( 0.215600);

        if( ATOMGROUP( LOCATION) == "C")
        {
          SetDisplacedSolventVolume( 16.44);
          SetRadius( 1.58);
          SetBoundHydrogen( 0);
        }
        else if( ATOMGROUP( LOCATION) == "CH")
        {
          SetDisplacedSolventVolume( 21.59);
          SetRadius( 1.73);
          SetBoundHydrogen( 1);
        }
        else if( ATOMGROUP( LOCATION) == "CH2")
        {
          SetDisplacedSolventVolume( 26.74);
          SetRadius( 1.85);
          SetBoundHydrogen( 2);
        }
        else
        {
          SetDisplacedSolventVolume( 31.89);
          SetRadius( 1.97);
          SetBoundHydrogen( 3);
        }
      }
      else if
      (
          ATOMGROUP( LOCATION) == "N" ||
          ATOMGROUP( LOCATION) == "NH" ||
          ATOMGROUP( LOCATION) == "NH2" ||
          ATOMGROUP( LOCATION) == "NH3"
      )
      {
        SetA1( 12.212600);
        SetA2( 3.132200);
        SetA3( 2.012500);
        SetA4( 1.166300);
        SetB1( 0.005700);
        SetB2( 9.893300);
        SetB3( 28.997500);
        SetB4( 0.582600);
        SetC( -11.529000);

        if( ATOMGROUP( LOCATION) == "N")
        {
          SetDisplacedSolventVolume( 2.49);
          SetRadius( 0.84);
          SetBoundHydrogen( 0);
        }
        else if( ATOMGROUP( LOCATION) == "NH")
        {
          SetDisplacedSolventVolume( 7.64);
          SetRadius( 1.22);
          SetBoundHydrogen( 1);
        }
        else if( ATOMGROUP( LOCATION) == "NH2")
        {
          SetDisplacedSolventVolume( 12.79);
          SetRadius( 1.45);
          SetBoundHydrogen( 2);
        }
        else
        {
          SetDisplacedSolventVolume( 17.94);
          SetRadius( 1.62);
          SetBoundHydrogen( 3);
        }
      }
      else if
      (
          ATOMGROUP( LOCATION) == "O" ||
          ATOMGROUP( LOCATION) == "OH"
      )
      {
        SetA1( 3.048500);
        SetA2( 2.286800);
        SetA3( 1.546300);
        SetA4(  0.867000);
        SetB1( 13.277100);
        SetB2(  5.701100);
        SetB3( 0.323900);
        SetB4( 32.908900);
        SetC(   0.250800);

        if( ATOMGROUP( LOCATION) == "O")
        {
          SetDisplacedSolventVolume( 9.13);
          SetRadius( 1.30);
          SetBoundHydrogen( 0);
        }
        else
        {
          SetDisplacedSolventVolume( 14.28);
          SetRadius( 1.5);
          SetBoundHydrogen( 1);
        }
      }
      else if
      (
         ATOMGROUP( LOCATION) == "S" ||
         ATOMGROUP( LOCATION) == "SH"
      )
      {
        SetA1( 6.905300);
        SetA2( 5.203400);
        SetA3( 1.437900);
        SetA4( 1.586300);
        SetB1( 1.467900);
        SetB2( 22.215100);
        SetB3( 0.253600);
        SetB4( 56.172000);
        SetC( 0.866900);

        if( ATOMGROUP( LOCATION) == "S")
        {
          SetDisplacedSolventVolume( 19.86);
          SetRadius( 1.68);
          SetBoundHydrogen( 0);
        }
        else
        {
          SetDisplacedSolventVolume( 25.10);
          SetRadius( 1.81);
          SetBoundHydrogen( 1);
        }
      }
      else if
      (
         ATOMGROUP( LOCATION) == "SE"
      )
      {
        SetA1( 17.000600);
        SetA2( 5.819600);
        SetA3( 3.973100);
        SetA4( 4.354300);
        SetB1( 2.409800);
        SetB2( 0.272600);
        SetB3( 15.237200);
        SetB4( 43.816300);
        SetC( 2.840900);
        SetDisplacedSolventVolume( 28.73);
        SetRadius( 1.90);
        SetBoundHydrogen( 0);
      }
      else
      {
        std::string grouptype( ATOMGROUP( LOCATION));
        BCL_Message( util::Message::e_Standard, " Something is wrong, the Atomtype is: " + util::Format()( grouptype));
      }
    }

    void FormFactorParameters::ShowValues()
    {
      BCL_Message( util::Message::e_Standard, " C1 scaling parameter : " + util::Format()(  m_ExcludedVolumeParameter ));
      BCL_Message( util::Message::e_Standard, " C2 scaling parameter : " + util::Format()(  m_HydrationShellParameter ));
      BCL_Message( util::Message::e_Standard, " Sasa parameter : " + util::Format()(  m_SolventAccessableSurfaceArea ));
      BCL_Message( util::Message::e_Standard, " Solvent parameter : " + util::Format()(  m_DisplacedSolventVolume ));
      BCL_Message( util::Message::e_Standard, " Radius : " + util::Format()(  m_Radius ));
      BCL_Message( util::Message::e_Standard, " Bound Hydrogen :http://www.sltrib.com " + util::Format()(  m_BoundHydrogen ));
      BCL_Message( util::Message::e_Standard, " A1 : " + util::Format()(  m_A1 ));
      BCL_Message( util::Message::e_Standard, " A2 : " + util::Format()(  m_A2 ));
      BCL_Message( util::Message::e_Standard, " A3 : " + util::Format()(  m_A3 ));
      BCL_Message( util::Message::e_Standard, " A4 : " + util::Format()(  m_A4 ));
      BCL_Message( util::Message::e_Standard, " B1 : " + util::Format()(  m_B1 ));
      BCL_Message( util::Message::e_Standard, " B2 : " + util::Format()(  m_B2 ));
      BCL_Message( util::Message::e_Standard, " B3 : " + util::Format()(  m_B3 ));
      BCL_Message( util::Message::e_Standard, " B4 : " + util::Format()(  m_B4 ));
      BCL_Message( util::Message::e_Standard, " C : " + util::Format()(  m_C ));
      BCL_Message( util::Message::e_Standard, " X : " + util::Format()(  m_X ));
      BCL_Message( util::Message::e_Standard, " Y : " + util::Format()(  m_Y ));
      BCL_Message( util::Message::e_Standard, " Z : " + util::Format()(  m_Z ));
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SaxsDebye::s_Instance
    (
      util::Enumerated< restraint::SasDebyeInterface>::AddInstance( new SaxsDebye())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SaxsDebye::SaxsDebye() :
      m_ShouldApproximateLoops( false),
      m_DetermineAnalyticNormFactor( false),
      m_ExcludedVolumeParameter( 1.0),
      m_HydrationShellParameter( 0.0),
      m_ShouldApproximateSideChains( true),
      m_ReducedExpData( util::ShPtr< storage::Vector< restraint::SasScatteringPoint> >()),
      m_Queue(),
      m_Program()
    {
    }

    //! @brief Constructor that takes a bool
    //! @param LOOPS bool value to represent loops that are not present in the protein model
    SaxsDebye::SaxsDebye
    (
      const CommandQueue &QUEUE,
      const bool LOOPS,
      const bool USE_REGULA_FALSI_APPROXIMATION,
      const float EXCLUDED_VOLUME_PARAMETER,
      const float HYDRATION_SHELL_PARAMETER,
      const bool SIDE_CHAIN_APPROXIMATION,
      const util::ShPtr< storage::Vector< restraint::SasScatteringPoint> > REDUCED_EXPERIMENTAL_DATA
    ) :
      m_ShouldApproximateLoops( LOOPS),
      m_DetermineAnalyticNormFactor( USE_REGULA_FALSI_APPROXIMATION),
      m_ExcludedVolumeParameter( EXCLUDED_VOLUME_PARAMETER),
      m_HydrationShellParameter( HYDRATION_SHELL_PARAMETER),
      m_ShouldApproximateSideChains( SIDE_CHAIN_APPROXIMATION),
      m_ReducedExpData( REDUCED_EXPERIMENTAL_DATA),
      m_Queue( QUEUE)
    {
      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_Saxs, util::CPPDataTypes::e_Float, m_Queue, std::string(), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

    //! @brief Clone function
    //! @return pointer to new SaxsDebye
    SaxsDebye *SaxsDebye::Clone() const
    {
      return new SaxsDebye( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SaxsDebye::GetAlias() const
    {
      static const std::string s_Name( "OpenCLSaxsDebye");
      return s_Name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SaxsDebye::GetClassIdentifier() const
    {
      // Get BCL Standardized Class name
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a proteinModel as input and calculates an intensity using the debye formula
    //! @param PROTEIN_MODEL
    //! @return the intensity for a given q value
    restraint::SasExperimentalAndCalculatedData SaxsDebye::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      util::ShPtr< restraint::SasScatteringData> sp_experimental_data;

      if( !this->GetReducedExperimentalData().IsDefined())
      {

        // get the experimental SAXS data
        sp_experimental_data = this->GetExperimentalData();

        // Verify you have experimental Saxs Data
        if( !sp_experimental_data.IsDefined())
        {
          // warn user and return empty data
          BCL_Message( util::Message::e_Critical, "No experimental SAXS data found, returning empty data");
          return restraint::SasExperimentalAndCalculatedData();
        }

      }
      else
      {

        restraint::SasScatteringData experimental_data;

        for
        (
          storage::Vector< restraint::SasScatteringPoint>::const_iterator
           exp_data_itr( this->GetReducedExperimentalData()->Begin()),
           exp_data_itr_end( this->GetReducedExperimentalData()->End());
          exp_data_itr != exp_data_itr_end;
          ++exp_data_itr
        )
        {
          experimental_data.PushBackScattering( *exp_data_itr);
        }

        // use the Clone to ShPtr to create a shared pointer to experimental data
        util::ShPtr< restraint::SasScatteringData> sp_reduced_data( util::CloneToShPtr( experimental_data));

        // set sp_experimental_data to the reduced data set
        sp_experimental_data = sp_reduced_data;
      }

      restraint::SasDebye saxs_debye_object
      (
        m_ShouldApproximateLoops,
        m_DetermineAnalyticNormFactor,
        m_ExcludedVolumeParameter,
        m_HydrationShellParameter,
        m_ShouldApproximateSideChains
      );

      saxs_debye_object.GetAtomsAndFormFactors( PROTEIN_MODEL);

       // create object to hold calculated data
      restraint::SasScatteringData calculated_data;

      calculated_data.AllocateScatteringMemory( sp_experimental_data->GetScatteringData().GetSize());

      const size_t number_of_atoms( saxs_debye_object.GetCoordinates().GetSize());

      // host_coords is pointer to a vector of 4 floats
      cl_float4 *host_coords;
      cl_float16 *host_params;

      // calloc allocates a block of memory for an array of num elements, each of them size bytes long
      // and initializes all its bits to zero

      // point host_coords to a block of memory the size of number of atoms by the size of cl_float4
      host_coords = ( cl_float4 *)calloc( number_of_atoms, sizeof( cl_float4));
      host_params = ( cl_float16 *)calloc( number_of_atoms, sizeof( cl_float16));

      FormFactorParameters data_set;

      // get the coordinates of each atom
      for( size_t row( 0), row_end( number_of_atoms); row < row_end; ++row)
      {
        host_coords[ row].s[0] = float( saxs_debye_object.GetCoordinates()( row)( 0));
        host_coords[ row].s[1] = float( saxs_debye_object.GetCoordinates()( row)( 1));
        host_coords[ row].s[2] = float( saxs_debye_object.GetCoordinates()( row)( 2));
        host_coords[ row].s[3] = float( 0);

        data_set.SetParameters
        (
          saxs_debye_object.GetAtomGroups(),
          saxs_debye_object.GetCoordinates(),
          saxs_debye_object.GetSASAPoint(),
          row,
          m_ExcludedVolumeParameter,
          m_HydrationShellParameter
        );

        host_params[ row] = data_set.GetParametersAsFloat16();
      }

      // initialize error number with success
      cl_int error_number = CL_SUCCESS;

      // bcl buffer - wrapper around opencl bindings
      // allocating memory on GPU, just an allocation
      Buffer device_coords( cl::Buffer( m_Queue.GetContext(), CL_FALSE, sizeof( cl_float4) * number_of_atoms, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

      // transfers host allocation of memory to the GPU
      error_number = m_Queue.enqueueWriteBuffer( device_coords, CL_FALSE, 0, sizeof( cl_float4) * number_of_atoms, host_coords);
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

      // allocating for parameters
      Buffer device_params( cl::Buffer( m_Queue.GetContext(), CL_FALSE, sizeof( cl_float16) * number_of_atoms, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

      // transfers host allocation of memory to the GPU
      error_number = m_Queue.enqueueWriteBuffer( device_params, CL_FALSE, 0, sizeof( cl_float16) * number_of_atoms, host_params);
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

      // get the number of q values
      const size_t number_q_values( sp_experimental_data->GetScatteringData().GetSize());

      // set up a vector of floats and doubles for the number of q values
      linal::Vector< float> q_values( number_q_values);
      linal::Vector< double> q_doubles( number_q_values);

      // reset row to 0
      size_t row( 0);

      // precalculate water factors and H for each q value
      // get the structure factor function for hydrogen

      // get the crommer mann constants for hydrogen
      static const math::SumFunction< restraint::SasDataParameters, double>
      h_form_factor( biol::GetAtomTypes().H->GetElementType()->GetStructureFactor(), 1.0, 0.0);

      // get the crommer mann constants for oxygen
      static const math::SumFunction< restraint::SasDataParameters, double>
      o_form_factor( biol::GetAtomTypes().O->GetElementType()->GetStructureFactor(), 1.0, 0.0);

      // initialize structure factors for water
      math::SumFunction< restraint::SasDataParameters, double> water_factors;
      math::SumFunction< restraint::SasDataParameters, double> h_factors;

      water_factors += double( 2.0) * h_form_factor;
      water_factors += o_form_factor;

      h_factors += h_form_factor;

      storage::Vector< float> water_factors_vector( number_q_values);

      storage::Vector< float> h_factors_vector( number_q_values);
      // iterate over experimental data to get q-values
      for
      (
        storage::Vector< restraint::SasScatteringPoint>::const_iterator
          data_itr( sp_experimental_data->GetScatteringData().Begin()),
          data_itr_end( sp_experimental_data->GetScatteringData().End());
        data_itr != data_itr_end;
        ++data_itr, ++row
      )
      {
        // variable to hold q-value both double and float forms
        q_values( row) = float( data_itr->GetQvalue());
        q_doubles( row) = data_itr->GetQvalue();

        restraint::SasDataParameters q_value( q_doubles( row));

        water_factors_vector( row) = float( water_factors( q_value));
        h_factors_vector( row) = float( h_factors( q_value));
      }

      //BCL_MessageDbg( " Water factor vector: " + util::Format()( water_factors_vector));
      //BCL_MessageDbg( " hydrogen factor vector: " + util::Format()( h_factors_vector));

      Vector< float> res_ff( number_of_atoms, m_Queue);

      // initialize
      Vector< float> inner_sum_matrix( number_of_atoms, m_Queue);

      // allocate memory to be transfered back
      linal::Vector< float> calculated_intensities( number_q_values);

      cl::Kernel inner_sum_kernel, inner_sum_kernel_zero;

      cl::Kernel res_ff_kernel( m_Program, "CalculateResFF", &error_number);

      BCL_Assert( error_number == CL_SUCCESS, "CalculateResFF arg error: " + opencl::Tools::ErrorString( error_number));
      if( GetTools().GetFirstCommandQueue().GetDevice( NULL).DeviceType( NULL) == CL_DEVICE_TYPE_GPU)
      {

        BCL_MessageDbg( " GPU Platform");
        // relies of implicit synchronization because of the warps
        inner_sum_kernel = cl::Kernel( m_Program, "InnerSumsOTFGPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

        inner_sum_kernel_zero = cl::Kernel( m_Program, "InnerSumsOTFZeroGPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));
      }
      else
      {
        BCL_MessageDbg( " CPU Platform");
        // does not rely on implicit synchronization
        inner_sum_kernel = cl::Kernel( m_Program, "InnerSumsOTFCPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

        inner_sum_kernel_zero = cl::Kernel( m_Program, "InnerSumsOTFZeroCPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));
      }

      // sum of intensities
      cl::Kernel intensity_sums_kernel( m_Program, "ReductionSumOfIntensity", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

      // rmsd group and worksize
      const cl_uint block_size = 128;
      const cl::NDRange block_dim( block_size);
      const cl::NDRange nr_atoms_worksize( Tools::RoundUp( block_size, number_of_atoms));
      const cl::NDRange offset;

      const size_t number_elements_partial_reduction( ( number_of_atoms / block_size) + 1);

      size_t q_index( 0);

      Vector< float> partial_reduction_output( number_elements_partial_reduction, m_Queue);

      if( q_values( q_index) == 0)
      {
        //BCL_MessageDbg( " q_values is zero ");

        error_number  = res_ff_kernel.setArg( 0, device_params);
        error_number  = res_ff_kernel.setArg( 1, q_values( q_index));
        error_number  = res_ff_kernel.setArg( 2, res_ff.GetData());
        error_number  = res_ff_kernel.setArg( 3, cl_uint( number_of_atoms));
        error_number  = res_ff_kernel.setArg( 4, h_factors_vector( q_index));
        error_number  = res_ff_kernel.setArg( 5, water_factors_vector( q_index));
        BCL_Assert( error_number == CL_SUCCESS, "ffs arg error: " + opencl::Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( res_ff_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        //BCL_MessageDbg( " gpu ffs: " + util::Format()( linal::Matrix< float>( res_ff.GetSize(), 1, res_ff.GetHostVector().Begin())));

        error_number  = inner_sum_kernel_zero.setArg( 0, device_coords);
        error_number |= inner_sum_kernel_zero.setArg( 1, res_ff.GetData());
        error_number |= inner_sum_kernel_zero.setArg( 2, q_values( q_index));
        error_number |= inner_sum_kernel_zero.setArg( 3, inner_sum_matrix.GetData());
        error_number |= inner_sum_kernel_zero.setArg( 4, cl_uint( number_of_atoms));
        BCL_Assert( error_number == CL_SUCCESS, "inner kernel zero arg error: " + opencl::Tools::ErrorString( error_number));

        // set
        error_number  = intensity_sums_kernel.setArg( 0, inner_sum_matrix.GetData());
        error_number |= intensity_sums_kernel.setArg( 1, cl_uint( number_of_atoms));
        error_number |= intensity_sums_kernel.setArg( 2, partial_reduction_output.GetData());
        error_number |= intensity_sums_kernel.setArg( 3, block_size * sizeof( float), 0);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

        // launch kernel
        error_number = m_Queue.enqueueNDRangeKernel( inner_sum_kernel_zero, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( intensity_sums_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        calculated_data.PushBackScattering( restraint::SasScatteringPoint( q_doubles( q_index), double( partial_reduction_output.GetHostVector().Sum()), 0.0));

        ++q_index;
      }

      for( ; q_index < number_q_values; ++q_index)
      {
        //BCL_MessageDbg( " q_values is not zero ");
        error_number  = res_ff_kernel.setArg( 0, device_params);
        error_number  = res_ff_kernel.setArg( 1, q_values( q_index));
        error_number  = res_ff_kernel.setArg( 2, res_ff.GetData());
        error_number  = res_ff_kernel.setArg( 3, cl_uint( number_of_atoms));
        error_number  = res_ff_kernel.setArg( 4, h_factors_vector( q_index));
        error_number  = res_ff_kernel.setArg( 5, water_factors_vector( q_index));
        BCL_Assert( error_number == CL_SUCCESS, "ffs arg error: " + opencl::Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( res_ff_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        //BCL_MessageDbg( " gpu ffs: " + util::Format()( linal::Matrix< float>( res_ff.GetSize(), 1, res_ff.GetHostVector().Begin())));

        // set
        error_number  = inner_sum_kernel.setArg( 0, device_coords);
        error_number |= inner_sum_kernel.setArg( 1, res_ff.GetData());
        error_number |= inner_sum_kernel.setArg( 2, q_values( q_index));
        error_number |= inner_sum_kernel.setArg( 3, inner_sum_matrix.GetData());
        error_number |= inner_sum_kernel.setArg( 4, cl_uint( number_of_atoms));
        BCL_Assert( error_number == CL_SUCCESS, "inner kernel arg error: " + opencl::Tools::ErrorString( error_number));

        // set
        error_number  = intensity_sums_kernel.setArg( 0, inner_sum_matrix.GetData());
        error_number |= intensity_sums_kernel.setArg( 1, cl_uint( number_of_atoms));
        error_number |= intensity_sums_kernel.setArg( 2, partial_reduction_output.GetData());
        error_number |= intensity_sums_kernel.setArg( 3, block_size * sizeof( float), 0);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

        // launch kernel
        error_number = m_Queue.enqueueNDRangeKernel( inner_sum_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( intensity_sums_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
        calculated_data.PushBackScattering( restraint::SasScatteringPoint( q_doubles( q_index), double( partial_reduction_output.GetHostVector().Sum()), 0.0));
      }

      restraint::SasExperimentalAndCalculatedData saxs_data( *sp_experimental_data, calculated_data);

      // Free the memory used by calloc
      delete [] host_coords;
      delete [] host_params;

      return saxs_data;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write errors out to
    bool SaxsDebye::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      // only return true if the opencl kernels are initialized; this might be false if -opencl Disable was passed
      return GetTools().HasCommandQueues();
    }

    //! @brief responsible for updating to a valid queue
    //! @param TOOLS opencl tools
    void SaxsDebye::UpdateQueue( Tools &TOOLS)
    {
      if( !TOOLS.HasCommandQueues())
      {
        return;
      }

      m_Queue = TOOLS.GetFirstCommandQueue();

      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_Saxs, util::CPPDataTypes::e_Float, m_Queue, std::string(), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SaxsDebye::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "performs saxs debye calculation using opencl for massively parallel architectures"
      );

      parameters.AddInitializer
      (
        "consider loops",
        "should loops be considered",
        io::Serialization::GetAgent( &m_ShouldApproximateLoops),
        "0"
      );
      parameters.AddInitializer
      (
        "analytic",
        "whether to determine the norm factor with regula falsi (1) or pythagorean approximation (0)",
        io::Serialization::GetAgent( &m_DetermineAnalyticNormFactor),
        "0"
      );
      parameters.AddInitializer
      (
        "excluded volume",
        "c1 adjustment parameter",
        io::Serialization::GetAgent( &m_ExcludedVolumeParameter),
        "1.0"
      );
      parameters.AddInitializer
      (
        "hydration shell",
        "c2 adjustment parameter",
        io::Serialization::GetAgent( &m_HydrationShellParameter),
        "0.0"
      );
      parameters.AddOptionalInitializer
      (
        "approximate_sidechains",
        "sum up form factor contribution on cb position",
        io::Serialization::GetAgent( &m_ShouldApproximateSideChains)
      );

      return parameters;
    }
  } // namespace opencl
} // namespace bcl
