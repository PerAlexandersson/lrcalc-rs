use lrcalc::schubert::{
    multiply_schubert, multiply_schubert_strings, strings_are_compatible, valid_permutation,
    LinearCombination, SchubertError,
};

fn main() {
    let mut args = std::env::args();
    let program = args.next().unwrap_or_else(|| "schubmult".to_string());
    let rest = args.collect::<Vec<_>>();
    match parse_args(&rest).and_then(run) {
        Ok(()) => {}
        Err(message) => {
            eprintln!("{program}: {message}");
            eprintln!("usage: {program} [-m] [-s] [-r rank] perm1 - perm2");
            std::process::exit(2);
        }
    }
}

struct Args {
    left: Vec<i32>,
    right: Vec<i32>,
    rank: i32,
    maple: bool,
    string_mode: bool,
}

fn parse_args(args: &[String]) -> Result<Args, String> {
    let mut maple = false;
    let mut string_mode = false;
    let mut rank = 0;
    let mut parts = Vec::new();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "-m" => {
                maple = true;
                index += 1;
            }
            "-s" => {
                string_mode = true;
                index += 1;
            }
            "-r" => {
                let value = args
                    .get(index + 1)
                    .ok_or_else(|| "missing value after -r".to_string())?;
                rank = value
                    .parse::<i32>()
                    .map_err(|_| format!("invalid rank value '{value}'"))?;
                if rank < 0 {
                    return Err("rank must be nonnegative".to_string());
                }
                index += 2;
            }
            token => {
                parts.push(token.to_string());
                index += 1;
            }
        }
    }
    if rank > 0 && string_mode {
        return Err("-s cannot be used with -r".to_string());
    }
    let [left, right] = parse_vector_pair(&parts)?;
    Ok(Args {
        left,
        right,
        rank,
        maple,
        string_mode,
    })
}

fn parse_vector_pair(args: &[String]) -> Result<[Vec<i32>; 2], String> {
    let mut parts = [Vec::new(), Vec::new()];
    let mut section = 0usize;
    for token in args {
        if token == "-" {
            section += 1;
            if section >= parts.len() {
                return Err("expected exactly one '-' separator".to_string());
            }
            continue;
        }
        for piece in token.split(',') {
            let piece = piece.trim_matches(|ch| matches!(ch, '(' | ')' | '[' | ']'));
            if piece.is_empty() {
                continue;
            }
            let value = piece
                .parse::<i32>()
                .map_err(|_| format!("invalid integer '{piece}'"))?;
            parts[section].push(value);
        }
    }
    if section != 1 {
        return Err("expected two vectors separated by '-'".to_string());
    }
    Ok(parts)
}

fn run(args: Args) -> Result<(), String> {
    let terms = if args.string_mode {
        if !strings_are_compatible(&args.left, &args.right) {
            return Err("incompatible strings".to_string());
        }
        multiply_schubert_strings(&args.left, &args.right).map_err(format_schubert_error)?
    } else {
        if !valid_permutation(&args.left) {
            return Err("perm1 not a valid permutation".to_string());
        }
        if !valid_permutation(&args.right) {
            return Err("perm2 not a valid permutation".to_string());
        }
        multiply_schubert(&args.left, &args.right, args.rank).map_err(format_schubert_error)?
    };

    if args.maple {
        print_maple_terms(&terms);
    } else {
        print_terms(&terms);
    }
    Ok(())
}

fn print_terms(terms: &LinearCombination) {
    for (key, coefficient) in terms {
        if *coefficient != 0 {
            println!("{coefficient}  ({})", format_vector(key));
        }
    }
}

fn print_maple_terms(terms: &LinearCombination) {
    print!("0");
    for (key, coefficient) in terms {
        if *coefficient == 0 {
            continue;
        }
        if *coefficient < 0 {
            print!("-{}*X[{}]", coefficient.unsigned_abs(), format_vector(key));
        } else {
            print!("+{coefficient}*X[{}]", format_vector(key));
        }
    }
    println!();
}

fn format_vector(values: &[i32]) -> String {
    values
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn format_schubert_error(error: SchubertError) -> String {
    match error {
        SchubertError::InvalidInput => "invalid input".to_string(),
        SchubertError::ArithmeticOverflow => "arithmetic overflow".to_string(),
    }
}
