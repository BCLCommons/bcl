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
#include "fold/bcl_fold_default_mutates.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "io/bcl_io_fixed_line_width_writer.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_default_mutates.cpp
  //!
  //! @author weinerbe
  //! @date Oct 10, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldDefaultMutates :
    public ExampleInterface
  {
  public:

    ExampleFoldDefaultMutates *Clone() const
    {
      return new ExampleFoldDefaultMutates( *this);
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
      // set flags for native
      fold::DefaultFlags::GetFlagNativeModel()->SetFlag();
      std::ostringstream error_stream;
      fold::DefaultFlags::GetFlagNativeModel()->GetParameterList().FirstElement()->SetParameter
      (
        AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"),
        error_stream
      );
      fold::DefaultFlags::GetFlagUseNativeSSEsAsPool()->SetFlag();

    ////////////////
    // operations //
    ////////////////

      // InitializeMutates
      BCL_MessageStd
      (
        "Initial mutate enums: " + util::Format()( fold::GetMutates().GetEnumCount())
      );
      fold::DefaultMutates::GetInstance().InitializeMutates();
      BCL_ExampleCheck( fold::GetMutates().GetEnumCount() != 0, true);
      BCL_MessageStd
      (
        "Final mutate enums: " + util::Format()( fold::GetMutates().GetEnumCount())
      );

      // ModifyMutateTree
      fold::MutateTree mutate_tree;
      fold::DefaultMutates::GetInstance().ModifyMutateTree( mutate_tree);
      BCL_ExampleCheck( mutate_tree.GetMutateProbabilities().IsEmpty(), false);
      BCL_MessageStd
      (
        "Added probabilities for " + util::Format()( mutate_tree.GetMutateProbabilities().GetSize()) + " mutate types"
      );

      // write the table
      std::stringstream string_stream;
      mutate_tree.CreateTable().WriteFormatted( string_stream);
      BCL_MessageVrb( string_stream.str());
      io::FixedLineWidthWriter flww;
      mutate_tree.WriteHelp( flww);
      BCL_Debug( flww.String());

      // unset the flags
      fold::DefaultFlags::GetFlagNativeModel()->UnsetFlag();
      fold::DefaultFlags::GetFlagUseNativeSSEsAsPool()->UnsetFlag();

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldDefaultMutates

  const ExampleClass::EnumType ExampleFoldDefaultMutates::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldDefaultMutates())
  );

} // namespace bcl
