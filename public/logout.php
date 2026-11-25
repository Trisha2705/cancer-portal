```php
<?php
session_start();

// Clear all session variables
$_SESSION = [];

// Delete the session cookie
if (ini_get('session.use_cookies')) {
    $params = session_get_cookie_params();

    setcookie(
        session_name(),
        '',
        time() - 42000,
        $params['path'],
        $params['domain'],
        $params['secure'],
        $params['httponly']
    );
}

// Destroy the session
session_destroy();
?>

<!DOCTYPE html>
<html lang="en">

<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <title>Signed out | CancerInsight</title>

    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"
          rel="stylesheet">

    <link rel="stylesheet" href="style/logout.css">
</head>

<body>

    <main id="logout-box">

        <div class="success-icon" aria-hidden="true">

            <svg viewBox="0 0 24 24"
                 width="32"
                 height="32"
                 fill="none"
                 stroke="currentColor"
                 stroke-width="2.5"
                 stroke-linecap="round"
                 stroke-linejoin="round">

                <polyline points="5 12.5 10 17.5 19 7.5"></polyline>

            </svg>

        </div>


        <div class="logout-header">

            <h1>You've been signed out</h1>

            <p class="subtitle">
                Thanks for using CancerInsight.
            </p>

        </div>


        <a href="login.php" class="return-button">
            Return to sign in
        </a>

    </main>

</body>
</html>
```