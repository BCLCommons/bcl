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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "chemistry/bcl_chemistry_small_molecule_qsar_storage_file.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_small_molecule_qsar_storage_file.cpp
  //!
  //! @author butkiem1
  //! @date August 1, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySmallMoleculeQsarStorageFile :
    public ExampleInterface
  {
  public:

    ExampleChemistrySmallMoleculeQsarStorageFile *Clone() const
    {
      return new ExampleChemistrySmallMoleculeQsarStorageFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////
    // construction //
    //////////////////

      // default constructor
      chemistry::SmallMoleculeQsarStorageFile smallmol_qsar_storage_default;

      // create a label that can initialize the object
      const std::string init_label
      (
        smallmol_qsar_storage_default.GetAlias()
        + "(filename=" + AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf") + ")"
      );

      util::Implementation< model::RetrieveDataSetBase> smallmol_qsar_storage( init_label);

    /////////////////
    // data access //
    /////////////////

      // load in a label, whose arguments the database will find
      io::IFStream input_label;
      BCL_ExampleMustOpenInputFile( input_label, AddExampleInputPathToFilename( e_Chemistry, "smallmolecule_qsar_storage_code.object"));
      util::ObjectDataLabel label( input_label);
      io::File::CloseClearFStream( input_label);

      BCL_MessageStd( "Read in label: " + label.ToString());

      // set the database up to find the feature and result, both the same in this toy example
      smallmol_qsar_storage_default.SelectFeatures( label);

      // check initializer
      BCL_ExampleIndirectCheck
      (
        util::Implementation< model::RetrieveDataSetBase>( init_label).IsDefined(),
        true,
        "Construction from label " + init_label
      );

      // class identifier
      BCL_ExampleCheck
      (
        smallmol_qsar_storage_default.GetClassIdentifier(),
        GetStaticClassName( smallmol_qsar_storage_default)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistrySmallMoleculeQsarStorageFile

  const ExampleClass::EnumType ExampleChemistrySmallMoleculeQsarStorageFile::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySmallMoleculeQsarStorageFile())
  );

} // namespace bcl
