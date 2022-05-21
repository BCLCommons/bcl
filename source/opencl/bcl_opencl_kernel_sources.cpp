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
#include "opencl/bcl_opencl_kernel_sources.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_file_in_search_path.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "opencl/bcl_opencl_kernel_source_file.h"
#include "opencl/bcl_opencl_tools.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

// path for opencl kernels
#if defined (__MINGW32__)
  #define BCL_KERNEL_PATH "opencl_kernels/"
#elif defined (__GNUC__)
  #define BCL_KERNEL_PATH "/dors/meilerlab/apps/bcl/opencl_kernels/rev_4941/"
#elif defined (_MSC_VER)
  #define BCL_KERNEL_PATH "../../../opencl_kernels/"
#endif

namespace bcl
{
  namespace opencl
  {

    const std::string &GetDefaultKernelsSourcePath()
    {
      static const std::string s_default_kernels_path
      (
        GetVersion().IsLicense() ?
          std::string( GetVersion().GetInstallPrefix() + "/opencl_kernels/")
          :
          std::string( BCL_KERNEL_PATH)
      );

      return s_default_kernels_path;
    }

    //! flag to change opencl kernels path
    util::ShPtr< command::FlagInterface> &KernelSources::GetKernelsSourcePathFlag()
    {
      static util::ShPtr< command::FlagInterface> s_kernels_path_flag
      (
        new command::FlagStatic
        (
          "opencl_kernels_path",
          "Path from which to retrieve opencl kernel (.cl) files",
          command::Parameter
          (
            "path",
            "relative or absolute path to a directory containing opencl kernels",
            command::ParameterCheckFileInSearchPath
            (
              "opencl_kernels",
              GetDefaultKernelsSourcePath(),
              io::Directory::e_Dir
            ),
            ""
          )
        )
      );

      return s_kernels_path_flag;
    }

    //! flag for path to to change opencl kernel binaries
    util::ShPtr< command::FlagInterface> &KernelSources::GetKernelsBinaryPathFlag()
    {
      static util::ShPtr< command::FlagInterface> s_kernels_path_flag
      (
        new command::FlagStatic
        (
          "opencl_kernels_binary_path",
          "if path is defined, all compiled binaries are stored and retrieved in the given directory after compilation",
          command::Parameter
          (
            "path",
            "relative or absolute path to read/write compiled binaries from/to",
            ""
          ),
          // This signal must be made after all opencl flags have been set
          // TODO: add SetSignal to flag interface so that this can be setup from where the flags are added, rather
          // than forcing it to be here
          &Tools::UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag
        )
      );

      return s_kernels_path_flag;
    }

    KernelSources::ProgramMap &KernelSources::GetBuildPrograms()
    {
      static ProgramMap s_build_programs;
      return s_build_programs;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    KernelSources::KernelSources() :
      e_EuclideanDistance            ( AddEnum( "EuclideanDistance" , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "euclidean_distance.cl"  )))),
      e_SvmKernelFunctions           ( AddEnum( "SvmKernelFunctions", util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "svm_kernel_functions.cl")))),
      e_InsertionSort                ( AddEnum( "InsertionSort"     , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "insertion_sort.cl"      )))),
      e_RMSD                         ( AddEnum( "RMSD"              , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "rmsd.cl"                )))),
      e_SequentialMinimalOptimization( AddEnum( "SMO"               , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "svr_kernels.cl"         )))),
      e_KappaNearestNeighbor         ( AddEnum( "KNN"               , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "knn_kernels.cl"         )))),
      e_ArtificialNeuralNetwork      ( AddEnum( "ANN"               , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "ann_kernels.cl"         )))),
      e_Linal                        ( AddEnum( "Linal"             , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "linal.cl"               )))),
      e_Quality                      ( AddEnum( "Quality"           , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "quality.cl"             )))),
      e_Saxs                         ( AddEnum( "SAXS"              , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "saxs_kernels.cl"        )))),
      e_Buffer                       ( AddEnum( "Buffer"            , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "buffer.cl"              )))),
      e_ClusterCenters               ( AddEnum( "ClusterCenters"    , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "cluster_centers.cl"     ))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KernelSources::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add the Opencl kernel path to filename
    //! @param FILE_NAME filename of kernel
    //! @return kernel_path/FILE_NAME
    std::string KernelSources::AddKernelPath( const std::string &FILE_NAME)
    {
      const std::string resolved_filename
      (
        util::GetRuntimeEnvironment().ResolveFileName( GetKernelsSourcePathFlag()->GetFirstParameter()->GetValue() + FILE_NAME)
      );

      BCL_Assert
      (
        !resolved_filename.empty(),
        "unable to resolve filename: " + GetKernelsSourcePathFlag()->GetFirstParameter()->GetValue() + FILE_NAME
      );

      return resolved_filename;
    }

    //! @brief compile a kernel and return a binary
    //! @param KERNEL the kernel to compile
    //! @param PRECISION the precision to use
    //! @param QUEUE the command queue to program will run in
    //! @param OPTIONS compiler options
    //! @param ERROR_PTR location to store the ERROR at
    //! @return the cl::Program, either the stored already compiled version, or the newly compiled version
    cl::Program
    KernelSources::Compile
    (
      const KernelSource &KERNEL,
      const util::CPPDataTypes::Types &PRECISION,
      const CommandQueue &QUEUE,
      const std::string &OPTIONS,
      cl_int *ERROR_PTR
    )
    {
      cl_int error_number( CL_SUCCESS);

      // program
      cl::Program program;

      // device the program will run on
      const Device device( QUEUE.GetDevice( &error_number));
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number);
        return program;
      }

      // platform
      const Platform platform( device.GetPlatform( &error_number));
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number);
        return program;
      }

      // context
      const Context context( QUEUE.GetContext( &error_number));
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number);
        return program;
      }

      // source identifier
      std::string identifier( ( *KERNEL)->GetIdentifier());
      // identifier refined with platform and device name
      if( !identifier.empty())
      {
        identifier = platform.Name() + device.Name() + OPTIONS + identifier;
        identifier += '_' + util::CPPDataTypes::GetCPPDatatypeName( PRECISION) + ".ptx";
        identifier = util::RemoveSpacesFromString( identifier);
      }

      // check if for given identifier, a program was already build
      ProgramMap::const_iterator itr_ident( GetBuildPrograms().find( identifier));
      if( itr_ident != GetBuildPrograms().end())
      {
        std::map< CommandQueue, cl::Program>::const_iterator queue_itr( itr_ident->second.find( QUEUE));
        if( queue_itr != itr_ident->second.end())
        {
          return queue_itr->second;
        }
      }

      io::DirectoryEntry bin_file;
      if( GetKernelsBinaryPathFlag()->GetFirstParameter()->GetWasSetInCommandLine() && !identifier.empty())
      {
        // file for that kernel, platform and device and options
        bin_file = io::DirectoryEntry( GetKernelsBinaryPathFlag()->GetFirstParameter()->GetValue() + identifier);

        // if such a file exists, compilation might not be necessary
        if( bin_file.DoesExist())
        {
          // read the file and make a program out of it
          cl::Program::Binaries binaries;
          io::IFStream read;
          io::File::MustOpenIFStream( read, bin_file.GetFullName(), std::ios::binary);
          const std::string file_content( ( std::istreambuf_iterator< char>( read)), std::istreambuf_iterator< char>());
          io::File::CloseClearFStream( read);
          binaries.push_back( std::make_pair( file_content.c_str(), file_content.size()));
          std::vector< cl_int> binary_status( 1, CL_SUCCESS); // number of devices - error for each device
          program = cl::Program( context, std::vector< cl::Device>( 1, device), binaries, &binary_status, &error_number);

          // was program creation successful
          if( error_number == CL_SUCCESS)
          {
            for( size_t i( 0); i < binary_status.size(); ++i)
            {
              if( binary_status[ i] != CL_SUCCESS)
              {
                error_number = binary_status[ i];
                BCL_MessageCrt( "unable to load binary for device: " + util::Format()( i) + " with error: " + Tools::ErrorString( binary_status[ i]));
                break;
              }
            }

            // no problem loading binaries
            if( error_number == CL_SUCCESS)
            {
              // compile
              error_number = program.build( std::vector< cl::Device>( 1, device), OPTIONS.c_str());
              if( error_number == CL_SUCCESS)
              {
                BCL_MessageVrb
                (
                  "use opencl binary file instead of sources: " + bin_file.GetFullName()
                );
                GetBuildPrograms()[ identifier][ QUEUE] = program;
                return program;
              }
            }
            // some error occurred
            error_number = CL_SUCCESS;
            program = cl::Program();
          }

          // file could not be used
          BCL_MessageCrt
          (
            "unable to use opencl binary file, recompiling for: " + bin_file.GetFullName()
          );
          identifier.clear();
        }
      }

      // actual source of the kernel
      const std::string ccc_source( ( *KERNEL)->GetSource( PRECISION, device.Extensions( NULL)));
      if( ccc_source.empty())
      {
        Tools::AssignError( ERROR_PTR, CL_INVALID_KERNEL_DEFINITION);
        return program;
      }

      // create the program
      cl::Program::Sources source;
      source.push_back( std::make_pair( ccc_source.c_str()      , ccc_source.length()));
      program = cl::Program( context, source, &error_number);
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number);
        return program;
      }

      // compile
      error_number = program.build( std::vector< cl::Device>( 1, device), OPTIONS.c_str());

      // error in compilation
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number); // assign the build error
        BCL_MessageCrt( "build program error: " + Tools::ErrorString( error_number));
        BCL_MessageCrt( "build program source:\n" + ccc_source);
        std::string build_info;
        error_number = program.getBuildInfo( device, CL_PROGRAM_BUILD_LOG, &build_info);
        if( error_number != CL_SUCCESS)
        {
          BCL_MessageCrt( "get build info error: " + Tools::ErrorString( error_number));
        }
        else
        {
          BCL_MessageCrt( "build log: " + build_info);
        }
        program = cl::Program();
        return program;
      }

      // successful compilation - store the binary if identifier is not empty
      if( GetKernelsBinaryPathFlag()->GetFirstParameter()->GetWasSetInCommandLine() && !identifier.empty() && !bin_file.DoesExist())
      {
        Tools::LogPtx( program, bin_file.GetFullName());
      }

      // cache the program build
      if( !identifier.empty())
      {
        GetBuildPrograms()[ identifier][ QUEUE] = program;
      }

      // end
      return program;
    }

    //! @brief construct on access function for all KernelSources
    //! @return reference to only instances of KernelSources
    KernelSources &GetKernelSources()
    {
      return KernelSources::GetEnums();
    }

  } // namespace opencl

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< opencl::KernelSourceInterface>, opencl::KernelSources>;

  } // namespace util
} // namespace bcl
