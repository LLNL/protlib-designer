import pytest
from click.testing import CliRunner
from scripts.run_protlib_designer import run_protlib_designer


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.parametrize(
    'args',
    [
        [
            'example_data/trastuzumab_spm.csv',
            '5',
            '--min-mut',
            '1',
            '--max-mut',
            '4',
            '--output-folder',
            'test_output',
        ],
        [
            'example_data/trastuzumab_spm.csv',
            '5',
            '--min-mut',
            '1',
            '--max-mut',
            '4',
            '--output-folder',
            'test_output',
            '--forbidden-aa',
            'C',
            '--max-arom-per-seq',
            '2',
            '--dissimilarity-tolerance',
            '0.0',
            '--interleave-mutant-order',
            'True',
            '--force-mutant-order-balance',
            'True',
            '--schedule',
            '1',
            '--schedule-param',
            '2,1',
            '--weighted-multi-objective',
            'True',
            '--data-normalization',
            'True',
        ],
        [
            'example_data/trastuzumab_spm.csv',
            '5',
            '--min-mut',
            '1',
            '--max-mut',
            '4',
            '--output-folder',
            'test_output',
            '--forbidden-aa',
            'C',
            '--max-arom-per-seq',
            '2',
            '--dissimilarity-tolerance',
            '0.5',
            '--interleave-mutant-order',
            'True',
            '--force-mutant-order-balance',
            'True',
            '--schedule',
            '2',
            '--schedule-param',
            '2,1',
            "--objective-constraints",
            '1',
            "--objective-constraints-param",
            '0.1',
            '--weighted-multi-objective',
            'True',
            '--data-normalization',
            'True',
        ],
    ],
)
def test_run_protlib_designer(runner, args, tmp_path):

    # Modify the output folder to be inside the temporary folder
    args[args.index('--output-folder') + 1] = str(tmp_path)

    result = runner.invoke(run_protlib_designer, args)
    assert result.exit_code == 0
