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
    let args = expand_short_options(args, &['r']);
    let mut maple = false;
    let mut string_mode = false;
    let mut rank = 0;
    let mut parts = Vec::new();
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--" => {
                parts.extend(args[index + 1..].iter().cloned());
                break;
            }
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

fn expand_short_options(args: &[String], value_options: &[char]) -> Vec<String> {
    let mut expanded = Vec::with_capacity(args.len());
    for arg in args {
        if arg == "-" || !arg.starts_with('-') || arg.starts_with("--") || arg.len() <= 2 {
            expanded.push(arg.clone());
            continue;
        }

        for (offset, option) in arg[1..].char_indices() {
            if !option.is_ascii_alphabetic() {
                expanded.push(arg.clone());
                break;
            }
            expanded.push(format!("-{option}"));
            if value_options.contains(&option) {
                let value_start = 1 + offset + option.len_utf8();
                if value_start < arg.len() {
                    expanded.push(arg[value_start..].to_string());
                }
                break;
            }
        }
    }
    expanded
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

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| (*value).to_string()).collect()
    }

    #[test]
    fn parser_accepts_clustered_short_options() {
        let parsed = parse_args(&args(&["-ms", "0", "1", "-", "1", "0"])).unwrap();
        assert!(parsed.maple);
        assert!(parsed.string_mode);
        assert_eq!(parsed.left, vec![0, 1]);
        assert_eq!(parsed.right, vec![1, 0]);

        let parsed = parse_args(&args(&["-r4", "2", "1", "-", "2", "1"])).unwrap();
        assert_eq!(parsed.rank, 4);

        let parsed = parse_args(&args(&["--", "1", "2", "-", "1", "2"])).unwrap();
        assert_eq!(parsed.left, vec![1, 2]);
        assert_eq!(parsed.right, vec![1, 2]);
    }
}
