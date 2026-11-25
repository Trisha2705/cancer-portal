
<?php
session_start();

require_once __DIR__ . '/../includes/db.php';

$error = '';

if ($_SERVER['REQUEST_METHOD'] === 'POST') {

    // Retrieve form values
    $email = trim($_POST['email'] ?? '');
    $password = $_POST['password'] ?? '';
    $role = $_POST['role'] ?? '';

    // Find the user based on their selected role
    if ($role === 'patient') {

        $stmt = $pdo->prepare(
            "SELECT *
             FROM users
             WHERE email = ?
             AND role = 'patient'"
        );

        $stmt->execute([$email]);

    } elseif ($role === 'oncologist' || $role === 'researcher') {

        $stmt = $pdo->prepare(
            "SELECT u.*
             FROM users u
             JOIN professionals p ON p.user_id = u.user_id
             WHERE u.email = ?
             AND u.role = 'professional'
             AND p.professional_type = ?"
        );

        $stmt->execute([$email, $role]);

    } else {
        $stmt = null;
    }

    $user = $stmt ? $stmt->fetch(PDO::FETCH_ASSOC) : false;

    // Check whether the user exists and their password is correct
    // The database column is password_hash
    // Check whether the user exists and their password is correct
    // The database column is password_hash
    if ($user && password_verify($password, $user['password_hash'])) {

        // Regenerate session ID for security
        session_regenerate_id(true);

        // Store user information in the session
        $_SESSION['user_id'] = $user['user_id'];
        $_SESSION['role'] = $role;

        // Redirect based on the verified role
        if ($role === 'patient') {

            header('Location: patient-dashboard.php');

        } elseif ($role === 'oncologist') {

            header('Location: oncologist/oncologist-dashboard.php');

        } elseif ($role === 'researcher') {

            header('Location: researcher-dashboard.php');
    }

    exit();

} else {
    $error = 'Invalid email, password or account type.';
}
}
?>

<!DOCTYPE html>
<html lang="en">

<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <title>Login | CancerInsight</title>

    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>

    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"
          rel="stylesheet">

    <link rel="stylesheet" href="style/login.css">
</head>

<body>

    <div id="left-column">

        <div class="brand">
            <div class="brand-mark" aria-hidden="true"></div>
            <span class="brand-name">CancerInsight</span>
        </div>

        <div class="hero-text">
            <h2>Understanding Cancer through Genetic Insights</h2>

            <p>
                Explore genetic mutations, investigate their associations
                with different cancer types and gain predictive insights
                through an integrated healthcare information platform.
            </p>
        </div>

    </div>

    <div id="right-column">

        <div id="login-box">

            <div class="login-header">
                <h1>Welcome to CancerInsight</h1>
                <p class="subtitle">Sign in to continue</p>
            </div>

            <div id="role-selector"
                 role="group"
                 aria-label="Select your role">

                <button type="button"
                        class="role-button active"
                        data-role="patient"
                        aria-pressed="true">
                    Patient
                </button>

                <button type="button"
                        class="role-button"
                        data-role="oncologist"
                        aria-pressed="false">
                    Oncologist
                </button>

                <button type="button"
                        class="role-button"
                        data-role="researcher"
                        aria-pressed="false">
                    Researcher
                </button>

            </div>

            <form id="login-form" method="POST" action="login.php">

                <input type="hidden"
                       name="role"
                       id="role"
                       value="patient">

                <div class="field">

                    <label for="email">Email</label>

                    <input type="email"
                           name="email"
                           id="email"
                           autocomplete="email"
                           placeholder="name@example.com"
                           required>

                </div>

                <div class="field">

                    <div class="label-row">

                        <label for="password">Password</label>

                        <a href="forgot-password.html"
                           class="small-link">
                            Forgot password?
                        </a>

                    </div>

                    <input type="password"
                           name="password"
                           id="password"
                           autocomplete="current-password"
                           required>

                </div>

                <div id="login-error"
                     class="error-message"
                     role="alert"
                     <?= $error ? '' : 'hidden' ?>>

                    <?= htmlspecialchars($error, ENT_QUOTES, 'UTF-8') ?>

                </div>

                <div class="remember-row">

                    <input type="checkbox"
                           name="remember"
                           id="remember">

                    <label for="remember">Remember details</label>

                </div>

                <input type="submit"
                       value="Sign in"
                       class="submit-button">

            </form>

            <p class="help-text">
                Need an account?
                <a href="mailto:admin@cancerinsight.example">
                    Contact your administrator
                </a>
            </p>

        </div>

    </div>

    <script>

        const roleButtons = document.querySelectorAll('.role-button');
        const roleInput = document.getElementById('role');

        roleButtons.forEach((button) => {

            button.addEventListener('click', () => {

                // Remove active state from all buttons
                roleButtons.forEach((b) => {

                    b.classList.remove('active');
                    b.setAttribute('aria-pressed', 'false');

                });

                // Activate selected role
                button.classList.add('active');
                button.setAttribute('aria-pressed', 'true');

                // Update hidden form input
                roleInput.value = button.dataset.role;

            });

        });

    </script>

</body>
</html>
