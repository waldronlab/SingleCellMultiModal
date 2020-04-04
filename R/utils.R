.getAnswer <- function(msg, allowed)
{
    if (interactive()) {
        repeat {
            cat(msg)
            answer <- readLines(n = 1)
            if (answer %in% allowed)
                break
        }
        tolower(answer)
    } else {
        "n"
    }
}

