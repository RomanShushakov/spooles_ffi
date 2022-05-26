use crate::{SymmetryFlag, PivotingFlag, solve_using_lu};


#[test]
fn test_solve_using_lu_symm()
{
    let n_row = 4;
    let a_symm = vec![
        (0, 0, 5.0), 
        (0, 1, -4.0),
        (0, 2, 1.0),
        (0, 3, 0.0),
        (1, 1, 6.0),
        (1, 2, -4.0),
        (1, 3, 1.0),
        (2, 2, 6.0),
        (2, 3, -4.0),
        (3, 3, 5.0),
    ];
    let b = vec![
        (0, 0.0),
        (1, 1.0),
        (2, 0.0),
        (3, 0.0),
    ];
    let symmetry_flag = SymmetryFlag::Symmetric;
    let pivoting_flag = PivotingFlag::Pivoting;
    let expected = vec![1.6, 2.6, 2.4, 1.4];

    let x = solve_using_lu(n_row, a_symm, b, symmetry_flag, pivoting_flag).unwrap();

    for i in 0..expected.len()
    {
        assert!((x[i] - expected[i]).abs() < 1e-12);
    }
}


#[test]
fn test_solve_using_lu_nonsymm()
{
    let n_row = 3;
    let a_symm = vec![
        (0, 0, 3.0),
        (0, 1, -0.1),
        (0, 2, -0.2),
        (1, 0, 0.1),
        (1, 1, 7.0),
        (1, 2, -0.3),
        (2, 0, 0.3),
        (2, 1, -0.2),
        (2, 2, 10.0),
    ];
    let b = vec![
        (0, 7.85),
        (1, -19.3),
        (2, 71.4),
    ];
    let symmetry_flag = SymmetryFlag::NonSymmetric;
    let pivoting_flag = PivotingFlag::Pivoting;
    let expected = vec![3.0, -2.5, 7.0];

    let x = solve_using_lu(n_row, a_symm, b, symmetry_flag, pivoting_flag).unwrap();

    for i in 0..expected.len()
    {
        assert!((x[i] - expected[i]).abs() < 1e-12);
    }
}
